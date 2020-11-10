deposition <- function(x,vDep,rn=NULL,Sensor=NULL,Source=NULL,vDepSpatial=NULL,ncores=1){#,fracDepInside=0,vDepInside=0,ncores=1){

	# argument names better!
	RNG <- RNGkind(kind="L'Ecuyer-CMRG")

	on.exit(
		{
			RNGkind(kind=RNG[1])
		}
	)

  # convert old versions 
  sx <- as.character(substitute(x))
	x <- copy(x)
	setDT(x)
	switchNames(x)
  if(is.null(attr(x, "Version"))){
		warning(paste0("Object '", sx[min(length(sx), 2)], "' has not yet been converted to version 4.2+"))
		convert(x)
  }
  
	ModelInput <- attr(x,"ModelInput")
	Catalogs <- attr(x,"Catalogs")
	Cat.Path <- attr(x,"CatPath")

	if(is.null(rn)){
		Selrn <- x[,unique(rn)]
	} else {
		Selrn <- rn
	}
	if(is.null(Sensor)){
		SelSensor <- x[,unique(Sensor)]
	} else {
		SelSensor <- Sensor
	}
	if(is.null(Source)){
		SelSource <- x[,unique(Source)]
	} else {
		SelSource <- Source
	}
	Run <- x[rn %in% Selrn & Sensor %chin% SelSensor & Source %chin% SelSource,]

	N <- nrow(Run)
	Run[, vd_index := 1:.N]

	if(is.character(vDep)){
		vDep <- Run[,vDep,with=FALSE][[1]]
	} else {
		vDep <- rep(vDep,N)[seq_len(N)]
	}

	# vDepSpatial
	if(vdSpat <- !is.null(vDepSpatial)){
		# check names:
		nms <- names(vDepSpatial[[1]])
		if(!inherits(vDepSpatial[[2]],"Sources")){
			stop("Second list entry of argument 'vDepSpatial' must be of class 'Sources'!")
		}
		if(!all(nms %in% unique(vDepSpatial[[2]][,1]))){
			stop(paste(nms[!(nms %in% unique(vDepSpatial[[2]][,1]))],collapse=", "),": area not defined!")
		}
		
		for(i in nms){
			if(is.character(vDepSpatial[[1]][[i]])){
				vDepSpatial[[1]][[i]] <- Run[,vDepSpatial[[1]][[i]],with=FALSE][[1]]
			} else {
				vDepSpatial[[1]][[i]] <- rep(vDepSpatial[[1]][[i]],N)[seq_len(N)]
			}
		}
	}

	### procSensors:
	pSens <- procSensors(ModelInput[["Sensors"]])

	SourceNames <- unique(ModelInput[["Sources"]][,1])
	# fracDepList <- vector(length(SourceNames),mode="list")
	# names(fracDepList) <- SourceNames
	# vDepList <- fracDepList
	# if(length(fracDepInside)==1){
	# 	fracDepList[] <- fracDepInside
	# } else {
	# 	if(!all(SourceNames%in%names(fracDepInside)))stop("Missing sources names in argument fracDepInside!\nSource Names are: \"",paste(SourceNames,collapse="\", \""),"\"")
	# 	fracDepList[] <- fracDepInside[SourceNames]
	# }
	# if(length(vDepInside)==1){
	# 	vDepList[] <- vDepInside
	# } else {
	# 	if(!all(SourceNames%in%names(vDepInside)))stop("Missing sources names in argument vDepInside!\nSource Names are: \"",paste(SourceNames,collapse="\", \""),"\"")
	# 	vDepList[] <- vDepInside[SourceNames]
	# }
	# vDepList <- as.data.frame(vDepList)

	setkey(Run,rn,Sensor)
	if(sfIsRunning()){ncores <- length(sfGetCluster())}

	if(parl <- (ncores!=1 & N>1)){


		if(!sfIsRunning()){
			on.exit(sfStop())
			if(ncores<1)ncores <- parallel::detectCores(TRUE,FALSE)
			sfInit(TRUE,ncores)
			invisible(sfClusterEval(data.table::setDTthreads(1)))
			setDTthreads(ncores)
		}

		cl <- sfGetCluster()


		Key <- Run[,paste0(.I,collapse=","),by=.(Sensor,rn)]
		nk <- nrow(Key)

		# sort by LineSensor/PointSensor:
		sort_nindex <- character(nk)
		if(length(pSens$"LineSensors") > 0){
			add <- rep("1",nk)
			add[Key[,Sensor %chin% pSens$"LineSensors"]] <- "0"
			sort_nindex <- paste0(sort_nindex,add)
		}
		
		# sort by height:
		Sheights <- sapply(Key[,Sensor],function(x,y){
			min(y[x == y[, "Sensor Name"], "Sensor Height (m)"])
		},y=pSens$"Calc.Sensors")
		add <- rank(Sheights, ties.method = "min")
		sort_nindex <- paste0(sort_nindex,sprintf(paste0("%0",floor(log(max(add),10))+1,"i"),add))

		Key <- Key[order(sort_nindex)]
		####

		ilist <- lapply(strsplit(Key[,V1],",",fixed=TRUE),as.numeric)
		RunList <- lapply(ilist,function(x,y)y[x,], as.data.frame(Run))

		#### calculate and reorder
		cat("Parallel computing deposition corrected C/E...\nThis will take a moment...\n")
		if(vdSpat){
			OutList <- snow::clusterApplyLB(cl,RunList,.calcDep_Wrapper,Catalogs,
				Cat.Path,ModelInput[["Sources"]],pSens$"Calc.Sensors",vDep,vDepSpatial, "spatial")
		} else {
			OutList <- snow::clusterApplyLB(cl,RunList,.calcDep_Wrapper,Catalogs,
				Cat.Path,ModelInput[["Sources"]],pSens$"Calc.Sensors",vDep,vDepSpatial)
		}

		# Out <- rbindlist(OutList)[order(vd_index)][,vd_index := NULL]
		Out <- rbindlist(OutList)

	} else {

		OutList <- vector(N,mode="list")
		if(vdSpat){
			for(i in seq_len(N)){
				cat(i,"/",N,":\n")
				OutList[[i]] <- .calcDep_Spatial(Run[i,],Catalogs,Cat.Path,ModelInput[["Sources"]],pSens$"Calc.Sensors",vDep,vDepSpatial)
			}
		} else {
			for(i in seq_len(N)){
				cat(i,"/",N,":\n")
				OutList[[i]] <- .calcDep(Run[i,],Catalogs,Cat.Path,ModelInput[["Sources"]],pSens$"Calc.Sensors",vDep,vDepSpatial)
			}			
		}
		Out <- rbindlist(OutList)[, vd_index := seq_len(N)]
	}

	setnames(Out,paste0(names(Out),"_Dep"))
	Out <- merge(Run, Out, by.x = "vd_index", by.y = "vd_index_Dep")
	Out[, vd_index := NULL]

	setattr(Out,"vDep",list(vDep=vDep,vDepSpatial=vDepSpatial))
	setattr(Out,"class",c("deposition",class(Out)))
	# rm(Run,ModelInput,Catalogs,pSens)
	# {xalt <- matrix(0,2,3)
	# xneu <- gc()
	# while(abs(xalt[2,3]-xneu[2,3])>0){xalt<-xneu;xneu <- gc()}
	# }
	return(Out)
}
