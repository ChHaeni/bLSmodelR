.calcOutput <- function(SncRun, InputList, C.Path, cl = NULL, show_progress = TRUE,
    variables = 'CE') {

	### procSensors:
	InputList$Sensors <- procSensors(InputList[["Sensors"]])

	SensorNames <- unique(InputList$Sensors$"Calc.Sensors"[, "Sensor Name"])
	SncRun[,Calc.Sensor:=""]

	SonicList <- vector(length(SensorNames),mode="list")
	splitSensor <- strsplit(SncRun[,Sensor],",",fixed=TRUE)
	for(i in seq_along(SensorNames)){
        # use exact matching
        sensor_points <- unlist(strsplit(InputList$Sensors$PS_list[[SensorNames[i]]], split = ',', fixed = TRUE))
        find_points <- lapply(splitSensor, function(x) x %in% sensor_points)
		ind <- which(unlist(lapply(find_points, any)))
		if(length(ind) > 0){
			SonicList[[i]] <- SncRun[ind,]
            SonicList[[i]][, Calc.Sensor := unlist(lapply(ind, function(x) paste(
                splitSensor[[x]][find_points[[x]]], collapse = ',')))]
			SonicList[[i]][, Sensor := SensorNames[i]]
		}
	}
    # check is.null (sensor names may be missing)
	SncRun <- rbindlist(SonicList[!sapply(SonicList, is.null)])
	setattr(SncRun,"class",c("SncExt","data.table","data.frame"))
	Key <- SncRun[,paste0(.I,collapse=","),by=.(Sensor,rn)]

	Slist <- attr(InputList[["Sources"]], "SourceList")	
	SrcsNames <- names(Slist)
	for(i in SrcsNames){
		Slist[[i]] <- cbind(Feld=i,Slist[[i]],stringsAsFactors=FALSE)
	}
	Srcs <- rbindlist(Slist)
	setnames(Srcs,c("Plot","x","y","pid"))
	setkey(Srcs,Plot,pid)
	rm(Slist)

	nk <- nrow(Key)

	SncRun[, seed := sample.int(1E9, .N)]
		
	if (!is.null(cl) && (nk > 1)) {

		# sort by LineSensor/PointSensor:
		sort_nindex <- character(nk)
		if(length(InputList$Sensors$"LineSensors")){
			add <- rep("1",nk)
			add[Key[,Sensor %chin% InputList$Sensors$"LineSensors"]] <- "0"
			sort_nindex <- paste0(sort_nindex,add)
		}
		
		# sort by height:
		Sheights <- sapply(Key[, Sensor], function(x, y){
			min(y[x == y[, "Sensor Name"], "Sensor Height (m)"])
		}, y = InputList$Sensors$"Calc.Sensors")
		add <- rank(Sheights,ties.method="min")
		sort_nindex <- paste0(sort_nindex,sprintf(paste0("%0",floor(log(max(add),10))+1,"i"),add))

		Key <- Key[order(sort_nindex)]
		####

        # prepare parallel input list
        SncList <- SncRun[, I(lapply(Key[, I(strsplit(V1, split = ','))], \(i) 
                .SD[as.numeric(i), ]))]

        # fix DTthreads
        old_nthreads <- data.table::setDTthreads(1L)

		cat("\n***********\n")
		cat("Parallel computing C/E ratios.\n")
		cat("\n\t-> This might take some time depending on the calculation load!!! <-\n\n")
		OutList <- .clusterApplyLB(cl, SncList, .calcCE, InputList, Srcs, C.Path, 
            progress = show_progress, variables = variables)

        # fix DTthreads
        data.table::setDTthreads(old_nthreads)
	
	} else {
		
		cat("\nCalculating CE ratios.\n")
		
		OutList <- vector(nk,mode="list")

		for(i in seq(nk)){
			ilist <- as.numeric(unlist(strsplit(Key[i,V1],",",fixed=TRUE)))

			cat("\n***********\n")
			cat("Sonic Row:",Key[i,rn],"\n")
			cat("Sensor:",Key[i,Sensor],"\n")
			cat("Count:",i,"/",nk,"\n")

			OutList[[i]] <- .calcCE(SncRun[ilist, ], InputList, Srcs, C.Path, 
                variables = variables)

		}
	}

	Out <- rbindlist(OutList)
	setattr(Out,"CalcSteps",SncRun)
    # add gc/memory attribute
    cpu_mem <- .gather_mem(OutList)
    setattr(Out, 'cpu_mem', cpu_mem)

	invisible(Out)
}
