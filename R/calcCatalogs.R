.calcCatalogs <- function(SncRun, InputList, C.Path, cl = NULL) {

	
	cindex <- SncRun[,which(Cat.calc)]

	if(!length(cindex)&InputList[["Model"]][["TDonly"]]){
		cat("All TD catalogs existing...\n")
	}

	#### Calculate TD catalogs and CE:
	if(length(cindex)){

		Calc.names <- grep("^Calc.",names(SncRun),value=TRUE)
		names(Orig.names) <- Orig.names <- gsub("^Calc.","",Calc.names)
		Orig.names[c("ZSens","Su_Ustar","Sv_Ustar")] <- c("SensorHeight","sUu","sVu")

		for(i in cindex){
			
			SnRun <- SncRun[i,]
			SnRun[,(Orig.names) := mget(Calc.names)]

			if(InputList[["Model"]][["TDread"]]){
				if(SnRun[,Cat.exists]){
					CatRead <- SnRun[,Cat.Name]
					if(InputList[["Model"]][["overwriteTD"]]){
						CatNameNew <- CatRead	
					} else {
						CatNameNew <- createCatName(SnRun,format(Sys.time(),"%Y%m%d%H%M%S"))					
					}					
				} else {
					CatNameNew <- paste0(gsub("-[0-9]{5,}_[0-9]{5,}_[0-9]{5,}$","_",SnRun[,Cat.Name]),format(Sys.time(),"%Y%m%d%H%M%S"))
				}
				# rename Catalog correctly
				ind <- SncRun[,which(Cat.Name==Cat.Name[i])]
				SncRun[ind,Cat.Name := CatNameNew]
			} else {
				CatNameNew <- createCatName(SnRun,format(Sys.time(),"%Y%m%d%H%M%S"))
				# rename Catalog correctly
				SncRun[i,Cat.Name:=CatNameNew]
			}

			cat("\n***********\n")
			cat("Calculate TDs for Catalog:",CatNameNew,"\n")
			cat("TD Interval",which(i==cindex),"/",length(cindex),"\n")
			cat("Sensor Height (z-d) =",SnRun[,SensorHeight],"m\n")

			if(SnRun[,Cat.exists]){
				Catalog <- readCatalog(paste(C.Path,CatRead,sep="/"))
				## korrigiere uvw0 U0:
				Cat.N0 <- attr(Catalog,"N0")
				Catalog <- initializeCatalog(SnRun,Catalog=Catalog)
				uvwind <- (Cat.N0+1):SnRun[,N0]
				cat("~~~~~\nCalculating remaining",sprintf("%d",SnRun[,N0] - Cat.N0),"Trajectories\n~~~~~\n\n")
			} else {
				uvwind <- SnRun[,(1:N0)]
				####### Initialize Touchdown Catalog:
				Catalog <- initializeCatalog(SnRun)
			}

			uvw <- uvw0(Catalog)
			cat("~~~~~~~~ initialized cov bias:\n")
			cat("bias cov(u0,v0) =",-cov(uvw[,"u0"],uvw[,"v0"]),"\n")
			cat("bias cov(v0,w0) =",-cov(uvw[,"v0"],uvw[,"w0"]),"\n")
			cat("bias cov(u0,w0) =",-SnRun[,Ustar^2] - cov(uvw[,"u0"],uvw[,"w0"]),"\n")
			cat("~~~~~~~~\n")
			# calculate TDs:
			cat("\nCalculating TDs...\n")
			if (!is.null(cl)) {
                pindex <- parallel::clusterSplit(cl, uvwind)
                # fix DTthreads
                old_nthreads <- data.table::setDTthreads(1L)
				pList <- parallel::clusterApply(cl, pindex, coreModelWrapper, uvw[, "u0"], uvw[, "v0"], uvw[, "w0"], SnRun)
                # fix DTthreads
                data.table::setDTthreads(old_nthreads)
				attCat <- attributes(Catalog)
				for(p in 1:length(pindex)){
					Catalog <- rbind(Catalog,list(as.integer(pindex[[p]][pList[[p]]$Traj_IDOut]),pList[[p]]$TimeOut,pList[[p]]$xOut,pList[[p]]$yOut,pList[[p]]$wTDOut))
				}
				for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref")))setattr(Catalog,ac,attCat[[ac]])
				rm(list=c("pList","pindex","attCat"))
			} else {
				ListOut <- coreModel(uvw[uvwind,"u0"],uvw[uvwind,"v0"],uvw[uvwind,"w0"],SnRun[,SensorHeight],SnRun[,Ustar],SnRun[,L],SnRun[,Zo],SnRun[,bw],SnRun[,sUu],SnRun[,sVu],SnRun[,kv],SnRun[,C0],SnRun[,alpha],SnRun[,MaxFetch])
				attCat <- attributes(Catalog)
				Catalog <- rbind(Catalog,list(as.integer(uvwind[ListOut$Traj_IDOut]),ListOut$TimeOut,ListOut$xOut,ListOut$yOut,ListOut$wTDOut))
				for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref")))setattr(Catalog,ac,attCat[[ac]])
				rm(list=c("ListOut","attCat"))
			}

			# wTDcutoff:
			if(Catalog[wTD<InputList[["Model"]][["wTDcutoff"]],.N])Catalog[wTD<InputList[["Model"]][["wTDcutoff"]],wTD:=InputList[["Model"]][["wTDcutoff"]]]
			
			cat("\ndone\n")

			# write catalog:
			if(InputList[["Model"]][["overwriteTD"]]&&SnRun[,Cat.exists]){
				cat("Replacing existing TD Catalog\n")
				file.remove(file=paste0(C.Path,"/",CatNameNew))
			} else {
				cat("Writing new TD Catalog\n")
			}
			writeCatalog(Catalog,paste0(C.Path,"/",CatNameNew))
			SncRun[i,Cat.exists:=TRUE]
		}
	}

	return(invisible(SncRun))	
}
