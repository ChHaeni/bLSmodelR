.calcDep_Wrapper <- function(input_file, spatial = FALSE) {
    # TODO: wrap try around?
    input_list <- qs::qread(input_file, strict = TRUE)
    # TODO: wrap try around?
    file.remove(input_file)
	setDT(input_list[['RunList']])
	setkey(input_list[['RunList']], rn, Sensor)
	out <- vector("list", input_list[['RunList']][,.N])
	if (spatial) {
		for (rl in seq_along(out)) {
			out[[rl]] <- .calcDep_Spatial(
                input_list[['RunList']][rl, ], input_list[['Catalogs']], input_list[['Cat.Path']], 
                input_list[['Sources']], input_list[['Sensors']], input_list[['vDep']], 
                input_list[['vDepSpatial']]
            )
		}
	} else {
		for (rl in seq_along(out)) {
			out[[rl]] <- .calcDep(
                input_list[['RunList']][rl, ], input_list[['Catalogs']], input_list[['Cat.Path']], 
                input_list[['Sources']], input_list[['Sensors']], input_list[['vDep']], 
                input_list[['vDepSpatial']]
            )
		}
	}
    # gather mem usage
    cpu_mem <- .gather_mem(out)
	out <- data.table::rbindlist(out)[, vd_index := input_list[['RunList']][, vd_index]]
    setattr(out, 'cpu_mem', cpu_mem)
    # TODO: wrap try around?
    qs::qsave(out, res_file <- sub('parbls_input_', 'parbls_result_', input_file), preset = 'uncompressed')
    # return result file path (and try error states?)
    return(res_file)
}

.calcDep_Spatial <- function(Run, Catalogs, C.Path, Sources, CSnsrs, vd, vdSpatial) {

    # record gc/mem
    .record_now(start = TRUE)

	# on.exit(browser())
	vd_index <- Run[, vd_index]
	# browser()	
	if(Run[,N_TD > 0]){
		
		# N_TD_tot <- Run[,N_TD]
		# N_TD_sum <- 0

        # merge Catalogs with row
		Row <- Catalogs[Run][order(as.numeric(gsub(".*[.]([0-9]*)$","\\1",PointSensor)))]
        # TODO: get vDep from Run?
		vdep <- vd[vd_index]

		n <- nrow(Row)
		indCats <- Row[,Cat.Name]
		uniqueCats <- unique(indCats)
		nCats <- length(uniqueCats)
		cat("Sensor", Row[1, Sensor], "/ Source", Row[1, Source],"\n")
		Steps <- unique(round(seq(1,n,length.out=10))) %w/o% n
		N0 <- Row[1,N0]
		cName <- ""
		cSeed <- -1
		Src <- Sources[Sources[,1] %chin% Row[1,Source],]

		Ci <- vector(mode="list",length=n)
		UVW <- vector(mode="list",length=nCats)
		names(UVW) <- uniqueCats

        # data.frame-like or other - i.e. does spatial name column exist?
        if (!is.null(names(vdSpatial[[1]])) && all(names(vdSpatial[[1]]) %in% vdSpatial[[2]][, 1])) {
            # get vDepSpatial column names
            nms_Spatial <- names(vdSpatial[[1]])
            # get corresponding vDep and Sources obj.
            vd_Spatial <- Run[, I(mget(vdSpatial[[1]]))]
        } else {
            # get vDepSpatial column names
            nms_Spatial <- unlist(strsplit(Row[, get(vdSpatial[[1]])], split = ','))
            # get corresponding vDep and Sources obj.
            vd_Spatial <- Run[, I(mget(nms_Spatial))]
        }
        # get corresponding Sources
		Src_Spatial <- lapply(nms_Spatial,function(x,y)y[y[,1] %in% x,],y=vdSpatial[[2]])
        # name them
		names(Src_Spatial) <- names(vd_Spatial) <- nms_Spatial
        # exclude if NA values in vDep
        nms_Spatial <- nms_Spatial[!is.na(unlist(vd_Spatial))]

		for(i in seq_len(n)){

			# if(N_TD_sum == N_TD_tot){
			# 	# browser()
			# 	break
			# }

			# progressbar
			if(i %in% Steps)cat(sprintf("\r[%s%s%s] %1.0f%%",paste(rep(">",round(2*(i-1)/n*10)),collapse=""),"|",paste(rep(".",max(0,19-round(2*(i-1)/n*10))),collapse=""),(i-1)/n*100))
			# {xalt <- matrix(0,2,3)
			# xneu <- gc()
			# while(abs(xalt[2,3]-xneu[2,3])>0){xalt<-xneu;xneu <- gc()}
			# }
			# if(cName!=Row[i,Cat.Name] | cSeed!=Row[i,Subset_seed]){
			if(cName!=Row[i,Cat.Name]){
				cName <- Row[i,Cat.Name]
				cSeed <- Row[i,seed]
				Ctlg <- readCatalog(paste0(C.Path,"/",cName))			
				initializeCatalog(Row[i, ], Catalog = Ctlg)
				uvw <- uvw0(Ctlg)
				Cat.N0 <- attr(Ctlg,"N0")
				rotateCatalog(Ctlg,Row[i,WD])
				if(Cat.N0>N0){
					env <- globalenv()
				    oseed <- env$.Random.seed
				    set.seed(cSeed,kind="L'Ecuyer-CMRG")			
					takeSub <- sample.int(Cat.N0,N0)
			        if (is.null(oseed)) {
			            rm(list = ".Random.seed", envir = env)
			        } else {
			            assign(".Random.seed", value = oseed, envir = env)
			        }
					attCat <- attributes(Ctlg)
					indexNew <- 1:N0
					names(indexNew) <- takeSub
					Ctlg <- Ctlg[Traj_ID %in% takeSub,]
					Ctlg[,":="(Traj_ID=indexNew[as.character(Traj_ID)])]
					for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref","uvw0")))setattr(Ctlg,ac,attCat[[ac]])
					uvw <- uvw[takeSub,]
				}
				UVW[[cName]] <- uvw
			}

			tag_inside(Ctlg,Src,CSnsrs[chmatch(Row[i,PointSensor], CSnsrs[, "Point Sensor Name"]),])
			setkey(Ctlg,Traj_ID)
			# browser()
			Ctlg[,tagInsideQ := td_inside]

			if(Ctlg[,any(tagInsideQ)]){

				# N_TD_sum <- N_TD_sum + Ctlg[,sum(tagInsideQ)]

				### spatially inhomogeneous vdep:
				Ctlg[,vDep := vdep]
				for(j in nms_Spatial){
					tag_inside(Ctlg,Src_Spatial[[j]],CSnsrs[chmatch(Row[i,PointSensor],CSnsrs[, "Point Sensor Name"]),])
					Ctlg[(td_inside),vDep := vd_Spatial[[j]]]
				}
                # TODO: allow for deposition inside source
                #   -> set vdep to 0 inside and remove subsetting below
                #   -> set vdep to 0 before spatial vdep such that it can be set to non
                #           zero through spatial...
				Ctlg[,":="(
					wTD2 = 2/wTD
					)][,dep:=1][(!tagInsideQ), dep := exp(-vDep*wTD2)
				]

				Ci[[i]] <- Ctlg[Traj_ID %in% Traj_ID[(tagInsideQ)],.(CE = sum(as.numeric(tagInsideQ)*cumprod(dep)*wTD2)),by=Traj_ID]
			}
		}
		cat(paste0("\r[",paste0(rep(">",20),collapse=""),"] 100%\n"))

		# rm(Catalogs,CSnsrs,Ctlg,Src,uvw,Run,Row)
		
        .record_now()

		# browser()

		if(n>1){

			# check nulls
			is_null <- sapply(Ci,is.null)
			# # erweitere UVW:
			# is_nUVW <- sapply(UVW,is.null)
			# UVW[is_nUVW] <- UVW[!is_nUVW][1]

			# weights:
			wts <- rep(2,n)
			wts[c(1,n)] <- 1
			rwts <- wts/sum(wts)

			Cmean <- rep(0,n)
			Cmean[!is_null] <- sapply(Ci[!is_null],function(x)x[,sum(CE)])/N0
			# individual CEs
			CEs <- mapply(function(i,x,y){
				out <- rep(0,N0)
				if(i){
					out
				} else {
					out[x[,Traj_ID]] <- x[,CE]
					out - y
				}
			},i=is_null,x=Ci,y=Cmean,SIMPLIFY=TRUE)

			# Us per Cat
			Us <- do.call(cbind,lapply(UVW,function(x){
				x[,"u0"] - mean(x[,"u0"])
			}))
			# Um per Cat
			Um <- do.call(cbind,lapply(UVW,function(x){
				x[,"u0"]
			}))
			# Vs per Cat
			Vs <- do.call(cbind,lapply(UVW,function(x){
				x[,"v0"]
			}))
			# Ws per Cat
			Ws <- do.call(cbind,lapply(UVW,function(x){
				x[,"w0"]
			}))
			rwm <- outer(rwts,rwts)/sqrt(N0)
			Out <- data.frame(
				CE=sum(Cmean*rwts),
				CE_se=sqrt(sum(cov(CEs)*rwm)),
				uCE=sum(colMeans(Us[,indCats[!is_null],drop=FALSE]*CEs[,!is_null,drop=FALSE])*rwts[!is_null]),			
				uCE_se=sqrt(sum(cov(Us[,indCats]*CEs)*rwm)),			
				vCE=sum(colMeans(Vs[,indCats[!is_null],drop=FALSE]*CEs[,!is_null,drop=FALSE])*rwts[!is_null]),
				vCE_se=sqrt(sum(cov(Vs[,indCats]*CEs)*rwm)),
				wCE=sum(colMeans(Ws[,indCats[!is_null],drop=FALSE]*CEs[,!is_null,drop=FALSE])*rwts[!is_null]),
				wCE_se=sqrt(sum(cov(Ws[,indCats]*CEs)*rwm)),
				UCE=sum(colMeans(Um[,indCats[!is_null],drop=FALSE]*t(t(CEs[,!is_null,drop=FALSE]) + Cmean[!is_null]))*rwts[!is_null])
			)

		} else {
			# browser()
			Ci_v <- rep(0,N0)
			Ci_v[Ci[[1]][,Traj_ID]] <- Ci[[1]][,CE]
			Out <- data.frame(
				CE=(CE <- mean(Ci_v)),
				CE_se=sd(Ci_v)/sqrt(N0),
				uCE=mean((UVW[[1]][,"u0"]-mean(UVW[[1]][,"u0"]))*(Ci_v - CE)),
				uCE_se=sd((UVW[[1]][,"u0"]-mean(UVW[[1]][,"u0"]))*(Ci_v - CE))/sqrt(N0),
				vCE=mean(UVW[[1]][,"v0"]*(Ci_v - CE)),
				vCE_se=sd(UVW[[1]][,"v0"]*(Ci_v - CE))/sqrt(N0),
				wCE=mean(UVW[[1]][,"w0"]*(Ci_v - CE)),
				wCE_se=sd(UVW[[1]][,"w0"]*(Ci_v - CE))/sqrt(N0),
				UCE=mean(UVW[[1]][,"u0"]*Ci_v)
			)		
		}
	} else {
		Out <- data.frame(
			CE=0,
			CE_se=NA_real_,
			uCE=0,
			uCE_se=NA_real_,
			vCE=0,
			vCE_se=NA_real_,
			wCE=0,
			wCE_se=NA_real_,
			UCE=0
		)		
	}

    # record gc/mem
    setattr(Out, 'cpu_mem', .record_now(reset = TRUE))

	return(Out)	
}

.calcDep <- function(Run, Catalogs, C.Path, Sources, CSnsrs, vd, vdSpatial) {


    # record gc/mem
    .record_now(start = TRUE)

	# on.exit(browser())

	# browser()

	vd_index <- Run[, vd_index]

	if(Run[,N_TD > 0]){
		
		# N_TD_tot <- Run[,N_TD]
		# N_TD_sum <- 0
		Row <- Catalogs[Run][order(as.numeric(gsub(".*[.]([0-9]*)$","\\1",PointSensor)))]
		vdep <- vd[vd_index]

		n <- nrow(Row)
		uniqueCats <- Row[,unique(Cat.Name)]
		indCats <- Row[,Cat.Name]
		nCats <- length(uniqueCats)
		cat("Sensor", Row[1, Sensor], "/ Source", Row[1, Source],"\n")
		Steps <- unique(round(seq(1,n,length.out=10))) %w/o% n
		N0 <- Row[1,N0]
		cName <- ""
		cSeed <- -1
		Src <- Sources[Sources[,1] %chin% Row[1,Source],]

		Ci <- vector(mode="list",length=n)
		UVW <- vector(mode="list",length=nCats)
		names(UVW) <- uniqueCats

		for(i in seq_len(n)){

			# if(N_TD_sum == N_TD_tot){
			# 	# browser()
			# 	break
			# }

			# progressbar
			if(i %in% Steps)cat(sprintf("\r[%s%s%s] %1.0f%%",paste(rep(">",round(2*(i-1)/n*10)),collapse=""),"|",paste(rep(".",max(0,19-round(2*(i-1)/n*10))),collapse=""),(i-1)/n*100))
			# {xalt <- matrix(0,2,3)
			# xneu <- gc()
			# while(abs(xalt[2,3]-xneu[2,3])>0){xalt<-xneu;xneu <- gc()}
			# }
			# if(cName!=Row[i,Cat.Name] | cSeed!=Row[i,Subset_seed]){
			if(cName!=Row[i,Cat.Name]){
				cName <- Row[i,Cat.Name]
				cSeed <- Row[i, seed]
				Ctlg <- readCatalog(paste0(C.Path,"/",cName))			
				initializeCatalog(Row[i, ], Catalog = Ctlg)
				uvw <- uvw0(Ctlg)
				Cat.N0 <- attr(Ctlg,"N0")
				rotateCatalog(Ctlg,Row[i,WD])
				if(Cat.N0>N0){
					env <- globalenv()
				    oseed <- env$.Random.seed
				    set.seed(cSeed,kind="L'Ecuyer-CMRG")			
					takeSub <- sample.int(Cat.N0,N0)
			        if (is.null(oseed)) {
			            rm(list = ".Random.seed", envir = env)
			        } else {
			            assign(".Random.seed", value = oseed, envir = env)
			        }
					attCat <- attributes(Ctlg)
					indexNew <- 1:N0
					names(indexNew) <- takeSub
					Ctlg <- Ctlg[Traj_ID %in% takeSub,]
					Ctlg[,":="(Traj_ID=indexNew[as.character(Traj_ID)])]
					for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref","uvw0")))setattr(Ctlg,ac,attCat[[ac]])
					uvw <- uvw[takeSub,]
				}
				UVW[[cName]] <- uvw
			}

			tag_inside(Ctlg, Src, CSnsrs[chmatch(Row[i, PointSensor], CSnsrs[, "Point Sensor Name"]), ])
			setkey(Ctlg,Traj_ID)
			
			if(Ctlg[,any(td_inside)]){	

				# N_TD_sum <- N_TD_sum + Ctlg[,sum(td_inside)]		
				
				Ctlg[,":="(
					wTD2 = 2/wTD
					)][,dep:=1][(!td_inside), dep := exp(-vdep*wTD2)
				]

				Ci[[i]] <- Ctlg[Traj_ID %in% Traj_ID[(td_inside)],.(CE = sum(as.numeric(td_inside)*cumprod(dep)*wTD2)),by=Traj_ID]

			}
		}
		cat(paste0("\r[",paste0(rep(">",20),collapse=""),"] 100%\n"))

		# rm(Catalogs,CSnsrs,Ctlg,Src,uvw,Run,Row)
		
        .record_now()

		# browser()

		if(n>1){


			# check nulls
			is_null <- sapply(Ci,is.null)
			# # erweitere UVW:
			# is_nUVW <- sapply(UVW,is.null)
			# UVW[is_nUVW] <- UVW[!is_nUVW][1]

			# weights:
			wts <- rep(2,n)
			wts[c(1,n)] <- 1
			rwts <- wts/sum(wts)

			Cmean <- rep(0,n)
			Cmean[!is_null] <- sapply(Ci[!is_null],function(x)x[,sum(CE)])/N0
			# individual CEs
			CEs <- mapply(function(i,x,y){
				out <- rep(0,N0)
				if(i){
					out
				} else {
					out[x[,Traj_ID]] <- x[,CE]
					out - y
				}
			},i=is_null,x=Ci,y=Cmean,SIMPLIFY=TRUE)

			# Us per Cat
			Us <- do.call(cbind,lapply(UVW,function(x){
				x[,"u0"] - mean(x[,"u0"])
			}))
			# Um per Cat
			Um <- do.call(cbind,lapply(UVW,function(x){
				x[,"u0"]
			}))
			# Vs per Cat
			Vs <- do.call(cbind,lapply(UVW,function(x){
				x[,"v0"]
			}))
			# Ws per Cat
			Ws <- do.call(cbind,lapply(UVW,function(x){
				x[,"w0"]
			}))
			rwm <- outer(rwts,rwts)/sqrt(N0)
			Out <- data.frame(
				CE=sum(Cmean*rwts),
				CE_se=sqrt(sum(cov(CEs)*rwm)),
				uCE=sum(colMeans(Us[,indCats[!is_null],drop=FALSE]*CEs[,!is_null,drop=FALSE])*rwts[!is_null]),			
				uCE_se=sqrt(sum(cov(Us[,indCats]*CEs)*rwm)),			
				vCE=sum(colMeans(Vs[,indCats[!is_null],drop=FALSE]*CEs[,!is_null,drop=FALSE])*rwts[!is_null]),
				vCE_se=sqrt(sum(cov(Vs[,indCats]*CEs)*rwm)),
				wCE=sum(colMeans(Ws[,indCats[!is_null],drop=FALSE]*CEs[,!is_null,drop=FALSE])*rwts[!is_null]),
				wCE_se=sqrt(sum(cov(Ws[,indCats]*CEs)*rwm)),
				UCE=sum(colMeans(Um[,indCats[!is_null],drop=FALSE]*t(t(CEs[,!is_null,drop=FALSE]) + Cmean[!is_null]))*rwts[!is_null])
			)
			
		} else {
			# browser()
			Ci_v <- rep(0,N0)
			Ci_v[Ci[[1]][,Traj_ID]] <- Ci[[1]][,CE]
			Out <- data.frame(
				CE=(CE <- mean(Ci_v)),
				CE_se=sd(Ci_v)/sqrt(N0),
				uCE=mean((UVW[[1]][,"u0"]-mean(UVW[[1]][,"u0"]))*(Ci_v - CE)),
				uCE_se=sd((UVW[[1]][,"u0"]-mean(UVW[[1]][,"u0"]))*(Ci_v - CE))/sqrt(N0),
				vCE=mean(UVW[[1]][,"v0"]*(Ci_v - CE)),
				vCE_se=sd(UVW[[1]][,"v0"]*(Ci_v - CE))/sqrt(N0),
				wCE=mean(UVW[[1]][,"w0"]*(Ci_v - CE)),
				wCE_se=sd(UVW[[1]][,"w0"]*(Ci_v - CE))/sqrt(N0),
				UCE=mean(UVW[[1]][,"u0"]*Ci_v)
			)		
		}
	} else {
		Out <- data.frame(
			CE=0,
			CE_se=NA_real_,
			uCE=0,
			uCE_se=NA_real_,
			vCE=0,
			vCE_se=NA_real_,
			wCE=0,
			wCE_se=NA_real_,
			UCE=0
		)		
	}

    # record gc/mem
    setattr(Out, 'cpu_mem', .record_now(reset = TRUE))

	return(Out)	
}
