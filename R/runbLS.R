runbLS <- function(ModelInput, Cat.Path = NULL, ncores = NULL, TDonly = NULL, asDT = TRUE, simpleNames = asDT) {

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

	if(!is.null(TDonly))ModelInput[["Model"]]$TDonly <- TDonly
	if(!is.null(ncores))ModelInput[["Model"]]$ncores <- ncores

	Model <- ModelInput[["Model"]]

	if (Model[["TDwrite"]] || Model[["TDread"]]) {
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
	
	if(!is.null(ModelInput[["Sources"]])){
		ModelInput[["Sources"]] <- procSources(ModelInput[["Sources"]])
	} else if (!Model[['TDonly']]) {
		stop("No Source Area specified!\n")
    }

	tempCats <- FALSE
	if(!Model[["TDwrite"]]){
		if(is.null(Cat.Path))Cat.Path <- getwd()
		Model[["overwriteTD"]] <- FALSE
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
	
	if (inherits(ncores, 'cluster')) {
        # get clusters
		cl <- ncores
        # get number of cores, but set ncores to 1 to keep running at end
        Model[['ncores']] <- length(cl)
		ncores <- 1
	} else if (ncores > 1) {
		on.exit(
			{
		    parallel::stopCluster(cl)
			}, add = TRUE
		)
        # set up PSOCK clusters
        cl <- parallel::makePSOCKcluster(ncores)
        # set data.table threads to ncores on master
        data.table::setDTthreads(ncores)
	} else if (ncores != 1) {
		stop("Number of cores must be greater or equal to 1!")
	} else {
        cl <- NULL
	}

    if (!is.null(cl)) {
        # setup RNG stream
        parallel::clusterSetRNGStream(cl, sample.int(1e9, 6, TRUE))
        # set wd
		gwd <- getwd()
        parallel::clusterExport(cl, c('gwd', 'C.Path'), envir = environment())
        parallel::clusterCall(cl, setwd, gwd)
        # set data.table threads to 1 on slaves
        parallel::clusterEvalQ(cl, data.table::setDTthreads(1L))
    }

	ModelInput[["Model"]] <- Model

	Intervals <- prepareIntervals(ModelInput, Cat.Path, TRUE, ncores = cl)
	
	.calcCatalogs(Intervals, ModelInput, Cat.Path, cl)

	if(!ModelInput[["Model"]][["TDonly"]]){
		
		Out <- .calcOutput(Intervals, ModelInput, Cat.Path, cl)
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
		if(!is.na(ncores) && ncores>1) parallel::stopCluster(cl)
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

