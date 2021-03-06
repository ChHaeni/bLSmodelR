.calcCE <- function(SubRun,InputList,Srcs,C.Path){


	# SubRun <- SncRun[ilist,]
	# zu beginn on.exit(fehlerangabe!?)
	setDT(SubRun)

	N0 <- SubRun[1,N0]
		
	AllSensorNames <- unlist(strsplit(SubRun[,Calc.Sensor],",",fixed=TRUE))
	lasn <- length(AllSensorNames)
	sind <- chmatch(AllSensorNames, InputList$Sensors$"Calc.Sensors"[, "Point Sensor Name"])
	SensorPositions <- as.matrix(InputList$Sensors$"Calc.Sensors"[sind, c("x-Coord (m)", "y-Coord (m)")])
	rownames(SensorPositions) <- InputList$Sensors$"Calc.Sensors"[sind, "Point Sensor Name"]
	Sheight <- range(InputList$Sensors$"Calc.Sensors"[sind, "Sensor Height (m)"])
	if(length(sind) > 1){
		Sheight <- paste(sprintf("%1.2f", Sheight), collapse = " to ")
	} else {
		Sheight <- sprintf("%1.3f", Sheight[1])
	}
	SourceNames <- unlist(strsplit(SubRun[1,Source],",",fixed=TRUE))	
	Scalc <- copy(Srcs[SourceNames])
	setkey(Scalc,Plot,pid)
	Srange <- Scalc[,rbind(
		cbind(x=max(x)-min(SensorPositions[,"x-Coord (m)"]),y=max(y)-min(SensorPositions[,"y-Coord (m)"]))
		,cbind(x=min(x)-max(SensorPositions[,"x-Coord (m)"]),y=min(y)-max(SensorPositions[,"y-Coord (m)"]))
		)]
	# initialize C list
	Ci_sub <- vector(mode="list",length=lasn)
	names(Ci_sub) <- AllSensorNames
	Ci <- lapply(seq_along(SourceNames),function(x)Ci_sub)
	names(Ci) <- SourceNames
	# initialize uvw list
	UVW <- vector(mode="list",length=nrow(SubRun))
	# initialize Output
	Out <- SubRun[rep(1,length(SourceNames)), ][, SensorHeight := as.character(SensorHeight)]
	Out[, ":="(
			Source=SourceNames,SourceArea=attr(InputList[["Sources"]],"SAreas")[SourceNames],
			SensorHeight = Sheight,
			CE=0,CE_se=NA_real_,CE_lo=NA_real_,CE_hi=NA_real_,
			uCE=0,uCE_se=NA_real_,uCE_lo=NA_real_,uCE_hi=NA_real_,
			vCE=0,vCE_se=NA_real_,vCE_lo=NA_real_,vCE_hi=NA_real_,
			wCE=0,wCE_se=NA_real_,wCE_lo=NA_real_,wCE_hi=NA_real_,
			N_TD=0,TD_Time_avg=NA_real_,TD_Time_max=NA_real_,Max_Dist=NA_real_,UCE=0,N_Sensor=0,Calc.Sensor=NULL,seed=NULL)]
	setkey(Out,Source)
	# Max_Dist = max fetch inside Source area

	for(Row in seq(nsr <- nrow(SubRun))){


		SensorNames <- unlist(strsplit(SubRun[Row,Calc.Sensor],",",fixed=TRUE))
		
		if(nsr>1)cat(paste0("\n> Sub-Count: ",Row," / ",nsr," (",length(SensorNames),"/",lasn," sub-sensors)\n"))

		# read Catalog:
		Catalog <- readCatalog(paste0(C.Path,"/",SubRun[Row,Cat.Name]))

		## korrigiere Catalog wTD uvw0 U0:
		initializeCatalog(SubRun[Row,],Catalog=Catalog)
		uvw <- uvw0(Catalog)
		# subset?
		if(attr(Catalog,"N0")>N0){
			env <- globalenv()
	    oseed <- env$.Random.seed
	    set.seed(SubRun[Row, seed],kind="L'Ecuyer-CMRG")			
			takeSub <- sort(sample.int(attr(Catalog,"N0"),N0))
	        if (is.null(oseed)) {
	            rm(list = ".Random.seed", envir = env)
	        } else {
	            assign(".Random.seed", value = oseed, envir = env)
	        }
			attCat <- attributes(Catalog)
			indexNew <- 1:N0
			names(indexNew) <- takeSub
			Catalog <- Catalog[Traj_ID %in% takeSub,]
			Catalog[,":="(Traj_ID=indexNew[as.character(Traj_ID)])]
			for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref","uvw0")))setattr(Catalog,ac,attCat[[ac]])
			uvw <- uvw[takeSub,]
			setattr(Catalog,"uvw0",uvw)
			HeaderSub <- paste(attr(Catalog,"header"),"\n*** Subset of ",N0," Trajectories ***\n\n",sep="")
			class(HeaderSub) <- c("TDhead","character")
			setattr(Catalog,"N0",N0)
			setattr(Catalog,"header",HeaderSub)
		}
		# save uvw
		if(is.null(UVW[[Row]])) UVW[[Row]] <- uvw

		rotateCatalog(Catalog,SubRun[Row,WD])
		cat(paste0("\nMatching catalog: ",SubRun[Row,Cat.Name],"\n*****\n",attr(Catalog,"header"),"*****\n\n"))
		
		cat("\nGet TD inside source areas:\n")
		tagNear(Catalog,Srange)
		Catalog[,inside0:=inside]
		Catalog[,rn:=.I]

		cat("\n~~ Sensor Height (z-d) >",SubRun[Row,SensorHeight],"m < ~~\n")
		if(Catalog[,any(inside)]){
			combs <- expand.grid(SourceNames,SensorNames,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
			SensorNumbers <- chmatch(combs[,2],AllSensorNames)
			nc <- NROW(combs)
			Steps <- unique(round(seq(1,nc,length.out=11))) %w/o% nc
			for(cmb in seq(nc)){
				
				# progressbar
				if(cmb %in% Steps)cat(sprintf("\r[%s%s%s] %1.0f%%",paste(rep(">",round(2*(cmb-1)/nc*10)),collapse=""),"|",paste(rep(".",max(0,19-round(2*(cmb-1)/nc*10))),collapse=""),(cmb-1)/nc*100))
				
				Source <- combs[cmb,1]
				Sensor <- combs[cmb,2]
				SourceAreaRelative <- copy(Scalc[Source])[,":="(x=x-SensorPositions[Sensor,1],y=y-SensorPositions[Sensor,2])]
				Catalog[,inside:=inside0]
				tagNear(Catalog,SourceAreaRelative)
				
				if(Catalog[,any(inside)]){
					# tag Inside Source
					Catalog[,inside1:=inside]
					TDinside <- SourceAreaRelative[,
					{
						tagNear(Catalog[,inside:=inside1],.(x=x,y=y))
						cbind(ID=Catalog[(inside),rn],pnt.in.poly(Catalog[(inside),cbind(x,y)],cbind(x,y)))
					},by=pid][,sum(pip),by=ID]

					if(any(TDinside[,as.logical(V1)])){
						setkey(Catalog,rn)
						Catalog[TDinside,inside:=as.logical(V1)]

						# calc CE
						Ci[[Source]][[Sensor]] <- Catalog[(inside),.(CE = sum(2/wTD)),by=Traj_ID]
						# Max_Dist etc.
						setkey(Catalog,Traj_ID)
						Cat <- Catalog[Catalog[(inside),.(minTime=min(Time)),by=Traj_ID]][Time>=minTime,]
						Out[Source,":="(
							Max_Dist = max(Max_Dist,rotateCatalog(Cat,SubRun[Row,WD],back=TRUE)[,-min(x)],na.rm=TRUE),
							N_TD = N_TD + Catalog[,sum(inside)],
							TD_Time_avg = sum(TD_Time_avg,Catalog[(inside),-mean(Time)],na.rm=TRUE),
							TD_Time_max = max(TD_Time_max,Catalog[(inside),-min(Time)],na.rm=TRUE),
							N_Sensor = N_Sensor + 1
							)]
					}
				}
			}
			cat(paste0("\r[",paste0(rep(">",20),collapse=""),"] 100%\n"))
		}
	}
	rm(Catalog)
	# {xalt <- matrix(0,2,3)
	# xneu <- gc()
	# while(abs(xalt[2,3]-xneu[2,3])>0){xalt<-xneu;xneu <- gc()}
	# }
	cat("\nCalculating Source contributions\n\n")

	# TD_Time_avg
	if(Out[,any(N_TD>0)])Out[N_TD>0,TD_Time_avg:=TD_Time_avg/N_Sensor]
	Out[,N_Sensor:=NULL]

	# weights:
	wts <- rep(2,lasn)
	wts[c(1,lasn)] <- 1
	rwts <- wts/sum(wts)
	if(lasn > 1){
		# rwts korrekt sortieren -> welche Sensoren sind max/min node? -> welche pos haben diese?
		Sname <- SubRun[1,Sensor]
		AllSensorOrder <- InputList$Sensors$"Calc.Sensors"[
			match(
				AllSensorNames, 
				InputList$Sensors$"Calc.Sensors"[, "Point Sensor Name"]), 
			"Node"]
		rwts[AllSensorOrder] <- rwts
	}

	# uvw key:
	uvw_key <- SubRun[,.(SubSensor = unlist(strsplit(Calc.Sensor,",",fixed=TRUE))),by=.(row=1:nrow(SubRun))]
	setkey(uvw_key,"SubSensor")

	# U Matrix:
	UVW_mean <- lapply(UVW,function(x)colMeans(x))
	U_matrix <- matrix(sapply(1:nrow(SubRun),function(x)UVW[[x]][,"u0"] - UVW_mean[[x]]["u0"]),nrow=N0)
	V_matrix <- matrix(sapply(1:nrow(SubRun),function(x)UVW[[x]][,"v0"] - UVW_mean[[x]]["v0"]),nrow=N0)
	W_matrix <- matrix(sapply(1:nrow(SubRun),function(x)UVW[[x]][,"w0"] - UVW_mean[[x]]["w0"]),nrow=N0)

	for(i in Out[N_TD>0,Source]){

		# check NULL
		not_null <- !sapply(Ci[[i]],is.null)
		
		# Mean I:
		CE_mean <- sum(sapply(Ci[[i]][not_null],function(x)x[,sum(CE)])*rwts[not_null])/N0
		UCE_mean <- sum(sapply(AllSensorNames[not_null],function(x,y){
			ind <- y[[x]][,Traj_ID]
			sum(y[[x]][,CE]*UVW[[uvw_key[x,row]]][ind,"u0"])
		},y=Ci[[i]])*rwts[not_null])/N0

		# initialize C Matrix
		c_matrix <- matrix(-CE_mean,nrow=N0,ncol=lasn)

		# fill C Matrix:
		for(j in which(not_null)){
			c_matrix[Ci[[i]][[j]][,Traj_ID],j] <- Ci[[i]][[j]][,CE - CE_mean]
		}

		# outer weights
		orwts <- outer(rwts,rwts)

		# CE SE
		CE_se_add <- sqrt(sum(cov(c_matrix)*orwts)/N0)
		
		# uCE + SE
		uvwCE <- c_matrix*U_matrix[,uvw_key[AllSensorNames,row]]
		uCE_add <- sum(colSums(uvwCE)*rwts)/N0
		uCE_se_add <- sqrt(sum(cov(uvwCE)*orwts)/N0)

		# vCE + SE
		uvwCE <- c_matrix*V_matrix[,uvw_key[AllSensorNames,row]]
		vCE_add <- sum(colSums(uvwCE)*rwts)/N0
		vCE_se_add <- sqrt(sum(cov(uvwCE)*orwts)/N0)

		# uCE + SE
		uvwCE <- c_matrix*W_matrix[,uvw_key[AllSensorNames,row]]
		wCE_add <- sum(colSums(uvwCE)*rwts)/N0
		wCE_se_add <- sqrt(sum(cov(uvwCE)*orwts)/N0)

		# write results
		Out[i,CE := CE_mean]
		Out[i,UCE := UCE_mean]
		Out[i,uCE := uCE_add]
		Out[i,vCE := vCE_add]
		Out[i,wCE := wCE_add]
		# SE:
		Out[i,CE_se := CE_se_add]
		Out[i,uCE_se := uCE_se_add]
		Out[i,vCE_se := vCE_se_add]
		Out[i,wCE_se := wCE_se_add]

	}
	# {xalt <- matrix(0,2,3)
	# xneu <- gc()
	# while(abs(xalt[2,3]-xneu[2,3])>0){xalt<-xneu;xneu <- gc()}
	# }
	# Upper/Lower
	qt_lo <- Out[, qt(0.025, N0 - 1)]
	qt_hi <- Out[, qt(0.975, N0 - 1)]
	Out[,":="(
		CE_lo = CE + qt_lo * CE_se,
		CE_hi = CE + qt_hi * CE_se,
		uCE_lo = uCE + qt_lo * uCE_se,
		uCE_hi = uCE + qt_hi * uCE_se,
		vCE_lo = vCE + qt_lo * vCE_se,
		vCE_hi = vCE + qt_hi * vCE_se,
		wCE_lo = wCE + qt_lo * wCE_se,
		wCE_hi = wCE + qt_hi * wCE_se
		)]

	return(Out)

}
