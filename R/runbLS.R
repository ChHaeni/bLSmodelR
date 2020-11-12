runbLS <- function(ModelInput,Cat.Path=NULL,ncores=NULL,writeCsv=FALSE,asDT=TRUE,simpleNames=asDT){

	cat("\n**********************************************\n")
	cat("\nLocal Date/Time: ",format(Start <- Sys.time(),format="%d-%m-%Y %H:%M:%S"),"\n")
	cat("\n**********************************************\n")

	RNG <- RNGkind(kind="L'Ecuyer-CMRG")

	on.exit(
		{
			RNGkind(kind=RNG[1])
			cat("\nLocal Date/Time: ",format(Sys.time(),format="%d-%m-%Y %H:%M:%S"),"\n")
			cat("\n***bLSmodelR Run aborted!***\n")
		}
	)

	if(!inherits(ModelInput,"InputList")){
		stop("Argument ModelInput must be of class \"InputList\"")
	} else {
		ModelInputOriginal <- ModelInput
	}

	Version <- "4.2+"

	C.Path <- Cat.Path
	Cat.Path <- Cat.Path[[1]]

	if(!is.null(ncores))ModelInput[["Model"]]$ncores <- ncores

	Model <- ModelInput[["Model"]]

	if(Model[["TD.write"]]|Model[["TD.read"]]){
		Cat.Path <- switch(class(Cat.Path),
			"character" = {
				if(!dir.exists(Cat.Path)){
					stop("Cat.Path folder does not exist:\n\t-> input: \"",Cat.Path,"\"")
				}else{
					Cat.Path
				}
			},
			{cat("Select TD catalog directory...\n")
			tk_choose.dir(default = getwd(), caption = "Select TD catalog directory")}
		)
		if(is.na(Cat.Path)){cat("\n- tcltk interface canceled -\n");return(invisible(NULL))}
	} 
	
	if((writeCsv==TRUE|class(writeCsv)=="character")&is.null(ModelInput[["Sources"]])){
		stop("No Source Area specified!\n")
	}
	TD.only <- is.null(ModelInput[["Sources"]])
	if(!TD.only){
		ModelInput[["Sources"]] <- procSources(ModelInput[["Sources"]])
	}

	switch(class(writeCsv),
		"logical" = {
			if(writeCsv){
				cat("Select Output File...\n")
				Out.File <- tclvalue(tkgetSaveFile(initialfile = paste0("bLS_Output_",format(Sys.time(),"%y%m%d_%H%M"),".csv"), title = "Save Output..."))
			} else {
				Out.File <- ""
			}
		},
		"character" = {
			Out.File <- writeCsv
			writeCsv <- TRUE
		},
		stop("Can not interpret writeCsv argument!\n")
		)


	tempCats <- FALSE
	if(!Model[["TD.write"]]){
		if(is.null(Cat.Path))Cat.Path <- getwd()
		Model[["TD.overwrite"]] <- FALSE
		dir_CatPath <- dir(Cat.Path,pattern="Cat|Header",full.names=TRUE)
		tempCats <- TRUE
		on.exit(
			{
				new_dir <- dir(Cat.Path,pattern="Cat|Header",full.names=TRUE)
				remove.me <- new_dir %w/o% dir_CatPath
				assign(paste0("RemoveTemporary",format(Sys.time(),"%Y%m%d_%H%M%S")),remove.me,envir=.GlobalEnv)
			},
			add=TRUE
		)
	}


	if(!is.list(C.Path))C.Path <- Cat.Path

	ncores <- Model[["ncores"]]
	
	if(sfIsRunning()){
		parl <- TRUE
		cl <- sfGetCluster()
		sfClusterSetupRNG(seed=sample.int(1E9,6,TRUE))
		hosts <- sapply(cl,"[","host")
		if(length(unique(hosts))==1&length(C.Path)==1){
			C.Path <- rep(list(C.Path),length(hosts))
			names(C.Path) <- hosts
		}
		for(i in 1:length(hosts)){
			eval(parse(text=paste0("clusterEvalQ(cl[i],C.Path <- \"",C.Path[[hosts[[i]]]],"\")")))
		}
		### check access:
		write.table("",file=paste0(Cat.Path,"/testNodes"))
		readTest <- try(sfClusterEval(read.table(paste0(C.Path,"/testNodes"))))
		file.remove(paste0(Cat.Path,"/testNodes"))
		if(inherits(readTest,"try-error")){
			stop("\n\tCat.Path must be accessible to all nodes!\n\tIn case you did not choose to save your Catalog Files,\n\tsupply a shared Catalog folder anyway.\n\tTemporary catalogs need to be saved there...\n")
		}
		ncores <- 1
	}
	if(ncores > 1){
		on.exit(
			{
			sfStop()
			},add=TRUE
		)
		sfInit(parallel=TRUE,cpus=ncores,type="SOCK")
		sfClusterSetupRNG(seed=sample.int(1E9,6,TRUE))
		gwd <- getwd()
		sfExport("gwd","C.Path")
		invisible(sfClusterCall(setwd,gwd))
		invisible(sfClusterEval(data.table::setDTthreads(1)))
		parl <- TRUE
		cl <- sfGetCluster()
	} else if(ncores != 1) {
		stop("Number of cores must be bigger or equal to 1!")
	} else {
		parl <- FALSE
	}
	
	data.table::setDTthreads(ncores)

	ModelInput[["Model"]] <- Model



	Intervals <- prepareIntervals(ModelInput,Cat.Path,TRUE)
	
	# Probleme auf NC mit untenstehender paralellisierung!
	# if(parl && Intervals[, floor(sum(Cat.calc) / ncores) > 2]){
	# 	# calculate TDs:
	# 	cat("\nCalculating TD Catalogs...\n")
	# 	lbList <- clusterApplyLB(cl, lapply(Intervals[,which(Cat.calc)], function(i) Intervals[i, ]), .calcCatalogs, InputList = ModelInput, C.Path = Cat.Path, parl = FALSE)
	# 	Intervals[(Cat.calc), names(Intervals) := rbindlist(lbList)]
	# } else {
	# 	.calcCatalogs(Intervals,ModelInput,Cat.Path,parl)
	# }
	
	.calcCatalogs(Intervals,ModelInput,Cat.Path,parl, TD.only)

	if(!TD.only){
		
		Out <- .calcOutput(Intervals,ModelInput,Cat.Path)
		Intervals <- attr(Out,"CalcSteps")
		Out[,":="(
			Cat.Name = NULL,
			Cat.exists = NULL,
			Cat.calc = NULL
			)]
		SAn <- match("SourceArea",names(Out))
		SHn <- match("SensorHeight",names(Out))
		sn <- match("UCE",names(Out))
		setcolorder(Out,c(1,13,14,SHn,SAn:sn,2:12,15:(SHn-1),(SHn+1):(SAn-1)))

		Out[, c("Sensor_Swustar", "Calc.ZSens", "Calc.Ustar", "Calc.L", "Calc.Zo",
			"Calc.Su_Ustar", "Calc.Sv_Ustar", "Calc.bw", "Calc.C0", "Calc.kv", 
			"Calc.A", "Calc.alpha", "Calc.MaxFetch", "Calc.Sensor_Swustar", 
			"Calc.N0") := NULL]

		setorder(Out,rn,Sensor,Source)
		setkey(Intervals,rn,Sensor)

		if(writeCsv){
			cat("Writing Output to File:\n",Out.File,"\n")
			terr <- try(write.table(Out,file=Out.File,sep=";",na="",row.names=FALSE,append=FALSE))
			if(inherits(terr,"try-error")){
				cat("\n > Output File could not be saved! <\n")
			}
		}

		Catalogs <- Intervals[,{
			ps <- unlist(strsplit(Calc.Sensor,","))
			.(
				PointSensor = ps, 
				Cat.Name = rep(Cat.Name, length(ps)),
				seed = rep(seed, length(ps))
				)
		},by=.(rn,Sensor,Calc.Sensor)][,Calc.Sensor:=NULL]
		setkey(Catalogs,rn,Sensor,PointSensor)

		# set attributes:
		setattr(Out,"CatPath",Cat.Path)
		setattr(Out,"Catalogs",Catalogs)
		setattr(Out,"ModelInput", ModelInputOriginal)
		setattr(Out, "Version", Version)
		# setattr(Out,"sessionInfo",sessionInfo())
		setattr(Out,"class",c("bLSresult","data.table","data.frame"))

		# sort by row names
		Out <- Out[order(match(rn, row.names(ModelInput$Interval)))]
		
		if(any(grepl("_add$",names(Out))))setnames(Out,gsub("_add$","",names(Out)))
	
		if(!asDT){
			ID <- Out[,paste(match(rn,sort.int(unique(rn))),match(Sensor,sort(unique(Sensor))),match(Source,sort(unique(Source))),sep=".")]
			setDF(Out)
			rownames(Out) <- ID
		}

		if(!simpleNames){
			Out <- switchNames(Out,simple=FALSE)
		}

		cat("Finished Run!\nLocal Date/Time: ",format(End <- Sys.time(),format="%d-%m-%Y %H:%M:%S"),"\nDuration of model run:",sprintf("%1.1f",dur <- End - Start),attr(dur,"units"),"\n****************************************************\n****************************************************\n")
		setattr(Out,"ModelRunTime",dur)

	} else {
		### prep output:
		pSens <- procSensors(ModelInput[["Sensors"]])
		SensorNames <- unique(pSens$"Calc.Sensors"[,1])
		Intervals[,Calc.Sensor:=""]
		SonicList <- vector(length(SensorNames),mode="list")
		splitSensor <- strsplit(Intervals[,Sensor],",",fixed=TRUE)
		for(i in seq_along(SensorNames)){
			SonicList[[i]] <- Intervals[grepl(paste0("(^|,)",SensorNames[i],"([.]+[0-9]+|)(,|$)"),Sensor),]
			if(nrow(SonicList[[i]])){
				ind <- grep(paste0("(^|,)",SensorNames[i],"([.]+[0-9]+|)(,|$)"),Intervals[,Sensor])
				SonicList[[i]][,Calc.Sensor:=sapply(splitSensor[ind],function(x)paste0(grep(paste0("^",SensorNames[i],"([.]+[0-9]+|)$"),x,value=TRUE),collapse=","))]				
				SonicList[[i]][,Sensor:=SensorNames[i]]
			}
		}
		Intervals <- Out <- rbindlist(SonicList)
		setattr(Intervals,"class",c("SncExt","data.table","data.frame"))
		setkey(Intervals,rn,Sensor)
		Catalogs <- Intervals[,{
			ps <- unlist(strsplit(Calc.Sensor,","))
			.(
				PointSensor = ps,
				Cat.Name = rep(Cat.Name, length(ps)),
				seed = rep(NA_integer_, length(ps))
				)
		},by=.(rn,Sensor,Calc.Sensor)][,Calc.Sensor:=NULL]
		setkey(Catalogs,rn,Sensor,PointSensor)

		setattr(Out,"CalcSteps",Intervals)
		Out[,":="(
			Cat.Name = NULL,
			Cat.exists = NULL,
			Cat.calc = NULL
			)]
		setcolorder(Out,c(1,13,14,2:12,15:ncol(Out)))
		setorder(Out,rn,Sensor,Source)

		Out[, c("Sensor_Swustar", "Calc.ZSens", "Calc.Ustar", "Calc.L", "Calc.Zo",
			"Calc.Su_Ustar", "Calc.Sv_Ustar", "Calc.bw", "Calc.C0", "Calc.kv", 
			"Calc.A", "Calc.alpha", "Calc.MaxFetch", "Calc.Sensor_Swustar", 
			"Calc.N0") := NULL]

		setattr(Out,"CatPath",Cat.Path)
		setattr(Out,"Catalogs",Catalogs)
		setattr(Out,"ModelInput", ModelInputOriginal)
		setattr(Out, "Version", Version)
		# setattr(Out,"sessionInfo",sessionInfo())
		setattr(Out,"class",c("bLSresult","data.table","data.frame"))


		# sort by row names
		Out <- Out[order(match(rn, row.names(ModelInput$Interval)))]
		
		if(any(grepl("_add$",names(Out))))setnames(Out,gsub("_add$","",names(Out)))
	
		if(!asDT){
			ID <- Out[,paste(match(rn,sort.int(unique(rn))),match(Sensor,sort(unique(Sensor))),sep=".")]
			setDF(Out)
			rownames(Out) <- ID
		}
		
		if(!simpleNames){
			Out <- switchNames(Out,simple=FALSE)
		}

		cat("Finished Run!\nLocal Date/Time: ",format(End <- Sys.time(),format="%d-%m-%Y %H:%M:%S"),"\nDuration of model run:",sprintf("%1.1f",dur <- End - Start),attr(dur,"units"),"\n****************************************************\n****************************************************\n")
		setattr(Out,"ModelRunTime",dur)
	}

	

	on.exit({
		RNGkind(kind=RNG[1])
		if(!is.na(ncores) && ncores>1)sfStop()
		if(tempCats){
			new_dir <- dir(Cat.Path,pattern="Cat|Header",full.names=TRUE)
			remove.me <- new_dir %w/o% dir_CatPath
			assign(paste0("RemoveTemporary",format(Sys.time(),"%Y%m%d_%H%M%S")),remove.me,envir=.GlobalEnv)
		}
		rebuildCatListFile(Cat.Path)
	})

	{xalt <- matrix(0,2,3)
	xneu <- gc()
	while(abs(xalt[2,3]-xneu[2,3])>0){xalt<-xneu;xneu <- gc()}
	}

	return(invisible(Out))
}

