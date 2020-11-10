genInputList <- function(...,Tol.Zero=FALSE){
	
	Interval<-NULL;Sensors<-NULL;Sources<-NULL;Model<-NULL;Tolerances<-NULL

	ArgList <- list(...)
	if(length(ArgList)){
		for(i in seq_along(ArgList)){
			switch(class(ArgList[[i]]) %w/o% c("data.frame","list"),
				"Interval"= Interval <- rbind(Interval,ArgList[[i]]),
				"Sensors"= Sensors <- join(ArgList[[i]],Sensors),
				"Sources"= Sources <- join(ArgList[[i]],Sources),
				"Tolerances"= if(!is.null(Tolerances)) stop("Please provide only one Tolerances data.frame!\n") else Tolerances <- ArgList[[i]],
				"Model"= if(!is.null(Model)) stop("Please provide only one Model list!\n") else Model <- ArgList[[i]]
				)			
		}
	}	
	if(is.null(Interval)){
		stop("Please provide a data.frame of class 'Interval'!")
	}
	if(is.null(Sensors)){
		stop("Please provide a data.frame of class 'Sensors'!")
	}
	if(is.null(Model)){
		Model <- genModel()
	}
	if(is.null(Tolerances)|Tol.Zero){
		Tolerances <- genTolerances(Tol.Zero=Tol.Zero)
	}
	if(!Model$TDonly&&is.null(Sources)){
		warning("Argument TDonly in 'Model' data.frame is set to TRUE since no Sources are defined!")
		Model$TDonly <- TRUE
	}

	# extend to all combinations?
	Snames_compact <- paste(Snames <- unique(Sensors[, "Sensor Name"]), collapse = ",")
	Interval[Interval[,"Sensor Names (sep = \",\")"] %in% "","Sensor Names (sep = \",\")"] <- Snames_compact
	Interval[Interval[,"Source Names (sep = \",\")"] %in% "","Source Names (sep = \",\")"] <- paste(unique(Sources[,1]),collapse=",")
	# check consistency of sensor and source names
	if(!all(unique(unlist(strsplit(Interval[,"Sensor Names (sep = \",\")"],","))) %in% unique(Snames)))
		stop("All entries in column \"Sensor Names (sep = \",\")\" need a matching entry in column \"Sensor Name\" of the Sensors data.frame!\n")
	if(!Model$TDonly&any(Interval[,"Source Names (sep = \",\")"]==""))stop("Interval column \"Source Names (sep = \",\")\" must be specified to calculate C/E values!\n")
	if(!Model$TDonly&!all(unique(unlist(strsplit(Interval[,"Source Names (sep = \",\")"],","))) %in% Sources[,1]))stop("All entries in column \"Source Names (sep = \",\")\" need a matching entry in column \"SourceArea Name\" of the Sources data.frame!\n")
	# check sensor heights:
	checkZo <- ""
	Nms <- character(0)
	for(i in seq(nrow(Interval))){
		ind <- match(
				unique(unlist(strsplit(Interval[i, "Sensor Names (sep = \",\")"], ","))),
				Sensors[, "Sensor Name"]
				, nomatch = 0) %w/o% 0
		if(!any(is.na(c(Interval[i, "Zo [m]"], Interval[i, "d [m]"]))) && 
			any(wHt <- Interval[i, "Zo [m]"] + Interval[i, "d [m]"] > Sensors[ind, "Sensor Height (m)"])){
			Nms <- unique(c(Nms, Sensors[ind[wHt], "Sensor Name"]))
		}
	}
	if(length(Nms)>1){
		brks <- ceiling(length(Nms)/10)
		Nms2 <- character(brks)
		for(j in seq(brks)){
			if(j==brks&&length(Nms)%%10){
				Nms2[j] <- paste(Nms[seq(length(Nms)%%10)+(j-1)*10],collapse=", ")
			} else {
				Nms2[j] <- paste(Nms[seq(10)+(j-1)*10],collapse=", ")
			}
		}
		checkZo <- paste0(checkZo,"\nSensors:\n\t",paste(Nms2,collapse="\n\t"),"\nare below the provided (Zo + d) for some rows in the specified Interval data.frame!")
	} else if(length(Nms==1)) {
		checkZo <- paste0(checkZo,"\nSensor ",Sensors[ind[wHt], "Sensor Name"]," is below the provided (Zo + d) for some rows in the specified Interval data.frame!")
	}
	if(checkZo!="")stop(paste0(checkZo,"\n"),call.=FALSE)
	
	return(structure(list(Interval=Interval,Model=Model,Sensors=Sensors,Sources=Sources,Tolerances=Tolerances),class=c("InputList","list")))
}
