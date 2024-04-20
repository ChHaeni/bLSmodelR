runbLS <- function(ModelInput, Cat.Path = NULL, ncores = NULL, TDonly = NULL, 
    asDT = TRUE, simpleNames = asDT, memory_limit = NULL, show_progress = TRUE,
    variables = 'CE', skipCrossCheck = FALSE) {

	cat("\n**********************************************\n")
    cat("                MAIN MODEL RUN\n")
	cat("\nLocal Date/Time: ",
        format(Start <- Sys.time(), format = "%d-%m-%Y %H:%M:%S"), 
        "\n")
	cat("\n**********************************************\n")

	RNG <- RNGkind(kind="L'Ecuyer-CMRG")

	on.exit(
		{
			RNGkind(kind=RNG[1])
			cat("\nLocal Date/Time: ",format(Sys.time(),format="%d-%m-%Y %H:%M:%S"),"\n")
			cat("\n***bLSmodelR Run aborted!***\n")
		}
	)

    # check variable argument
    if (!all(variables %in% c('CE', 'wCE', 'uCE', 'vCE'))) {
        stop('argument "variables" should be a vector with any combination of',
            ' "CE", "uCE", "vCE" and "wCE"')
    }

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
        if (!is.character(Cat.Path)) {
            stop("Cat.Path argument is not of type character")
        } else if (!dir.exists(Cat.Path)) {
            stop("Cat.Path folder does not exist:\n\t-> input: \"", Cat.Path, "\"")
        }
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
		ncores <- length(cl)
        if (!is.null(memory_limit)) {
            warning('argument "memory_limit" has been provided,\nbut memory limitation',
                ' is not possible on already running workers.')
        }
	} else if (ncores > 1) {
		on.exit(
			{
		    if (exists('cl')) parallel::stopCluster(cl)
			}, add = TRUE
		)
        # set up PSOCK clusters
        cl <- .makePSOCKcluster(ncores, memory_limit = memory_limit)
        # set data.table threads to ncores on master
        # data.table::setDTthreads(ncores)
	} else if (ncores != 1) {
		stop("Number of cores must be greater or equal to 1!")
	} else {
        cl <- NULL
        if (!is.null(memory_limit)) {
            warning('argument "memory_limit" has been provided,\nbut memory limitation',
                ' is only possible when setting up workers for parallel processing.')
        }
	}

    if (!is.null(cl)) {
        cat('\n*** Parallel computation on', length(cl), 'cores ***\n\n')
        # setup RNG stream
        parallel::clusterSetRNGStream(cl, sample.int(1e9, 6, TRUE))
        # set wd
		gwd <- getwd()
        parallel::clusterExport(cl, 'gwd', envir = environment())
        parallel::clusterCall(cl, setwd, gwd)
        if (.is_recording()) {
            parallel::clusterEvalQ(cl, bLSmodelR:::.start_recording())
        }
        # set data.table threads to ncores on main
        data.table::setDTthreads(ncores)
        # set data.table threads to 1 on slaves
        parallel::clusterEvalQ(cl, data.table::setDTthreads(1L))
    }

	ModelInput[["Model"]] <- Model

	Intervals <- prepareIntervals(ModelInput, Cat.Path, TRUE, ncores = cl, 
        throttle = getOption('bls.throttle', 100), skipCrossCheck = skipCrossCheck)
	
	.calcCatalogs(Intervals, ModelInput, Cat.Path, cl)

	if(!ModelInput[["Model"]][["TDonly"]]){
		
		Out <- .calcOutput(Intervals, ModelInput, Cat.Path, cl, 
            show_progress = show_progress, variables = variables)
        cpu_mem <- attr(Out, 'cpu_mem')
		Intervals <- attr(Out,"CalcSteps")
		Out[,":="(
			Cat.Name = NULL,
			Cat.exists = NULL,
			Cat.calc = NULL
			)]
		SAn <- match("SourceArea", names(Out))
		SHn <- match("SensorHeight", names(Out))
		sn <- match("UCE", names(Out))
        kn1 <- match(c('kv', 'A', 'alpha'), names(Out))
        kn2 <- match(c('bw', 'C0'), names(Out))
        firstn <- c(1, 13, 14, SHn, SAn:sn, 2:12, kn1, kn2) 
        residn <- seq_along(names(Out)) %w/o% firstn
		setcolorder(Out, c(firstn, residn))

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
		
        # remove passthrough suffix
		setnames(Out, gsub("_add$", "", names(Out)))
	
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
        setattr(Out, 'cpu_mem', cpu_mem)

	} else {
		### prep output:
		pSens <- procSensors(ModelInput[["Sensors"]])
		SensorNames <- unique(pSens$"Calc.Sensors"[, 1])
		Intervals[, Calc.Sensor := ""]
        # get seeds from catalogs
        Intervals[, seed := {
            attr(
                readCatalog(file.path(Cat.Path, .BY[['Cat.Name']])),
                'UVWseed'
            )
        }, by = Cat.Name]

		SonicList <- vector(length(SensorNames), mode = "list")
		splitSensor <- strsplit(Intervals[, Sensor], ",", fixed = TRUE)
        for(i in seq_along(SensorNames)){
            # use exact matching
            sensor_points <- unlist(strsplit(pSens$PS_list[[SensorNames[i]]], split = ',', fixed = TRUE))
            find_points <- lapply(splitSensor, function(x) x %in% sensor_points)
            ind <- which(unlist(lapply(find_points, any)))
            if(length(ind) > 0){
                SonicList[[i]] <- Intervals[ind,]
                SonicList[[i]][, Calc.Sensor := unlist(lapply(ind, function(x) paste(
                    splitSensor[[x]][find_points[[x]]], collapse = ',')))]
                SonicList[[i]][, Sensor := SensorNames[i]]
            }
        }
		Intervals <- Out <- rbindlist(SonicList)
		setattr(Intervals,"class",c("SncExt","data.table","data.frame"))
		setkey(Intervals,rn,Sensor)

		Catalogs <- Intervals[, {
			ps <- unlist(strsplit(Calc.Sensor, ","))
			.(
				PointSensor = ps, 
				Cat.Name = rep(Cat.Name, length(ps)), 
				seed = rep(seed, length(ps))
				)
		}, by = .(rn, Sensor, Calc.Sensor)][, Calc.Sensor := NULL]
		setkey(Catalogs,rn,Sensor,PointSensor)

		setattr(Out,"CalcSteps",Intervals)
		Out[,":="(
			Cat.Name = NULL,
			Cat.exists = NULL,
			Cat.calc = NULL
			)]
		setcolorder(Out,c(1,13,14,2:12,15:ncol(Out)))
		setorder(Out,rn,Sensor,Source)

		setattr(Out,"CatPath",Cat.Path)
		setattr(Out,"Catalogs",Catalogs)
		setattr(Out,"ModelInput", ModelInputOriginal)
		setattr(Out, "Version", Version)
		# setattr(Out,"sessionInfo",sessionInfo())
		setattr(Out,"class",c("bLSresult","data.table","data.frame"))


		# sort by row names
		Out <- Out[order(match(rn, row.names(ModelInput$Interval)))]
		
		setnames(Out, gsub("_add$", "", names(Out)))
	
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

