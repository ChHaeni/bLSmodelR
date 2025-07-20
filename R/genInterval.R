genInterval <- function(Data=NULL,Ustar=0.25,L=-2000,Zo=0.01,sUu=2.5,sVu=2,sWu=1.25,z_sWu=2,WD=0,d=0,N0=5E4,MaxFetch=500,SensorNames="",SourceNames=""){
	sNames <- c(
		"Ustar","L","Zo","sUu","sVu","sWu","z_sWu","WD","d","N0","MaxFetch","SensorNames","SourceNames",
		"rn","SensorHeight","SourceArea","CE","CE_se","CE_lo","CE_hi","uCE","uCE_se",
		"uCE_lo","uCE_hi","vCE","vCE_se","vCE_lo","vCE_hi","wCE","wCE_se","wCE_lo","wCE_hi",
		"N_TD","TD_Time_avg","TD_Time_max","Max_Dist","UCE","Subset_seed","kv","A","alpha",
		"bw","C0","Sensor_Swustar","Sensor","Source"
		) 
	cNames <- c(
		"Ustar [m/s]","L [m]","Zo [m]","SigmaU/Ustar [-]","SigmaV/Ustar [-]","SigmaW/Ustar [-]",
		"SigmaW/Ustar Height [m]","WD [deg N]","d [m]","N0","MaxFetch","Sensor Names (sep = \",\")",
		"Source Names (sep = \",\")","Original Interval Row","Sensor Height [m above d]",
		"Source Area [m]","C/E","C/E SE","C/E Lower-95%CI","C/E Upper-95%CI","u'C'/E","u'C'/E SE",
		"u'C'/E Lower-95%CI","u'C'/E Upper-95%CI","v'C'/E","v'C'/E SE","v'C'/E Lower-95%CI",
		"v'C'/E Upper-95%CI","w'C'/E","w'C'/E SE","w'C'/E Lower-95%CI","w'C'/E Upper-95%CI",
		"Number of TD inside","TD Time avg [sec]","TD Time max [sec]","Maximum Distance [m]","uC/E",
		"Catalog Subset Seed","SigmaW/Ustar @ Sensor Height"
		)	
	sNames13 <- sNames[1:13]
	cNames13 <- cNames[1:13]
	
	if(!is.null(Data)){
		
		# get call
		mc <- as.list(match.call()[-1])
		mc <- mc[!(names(mc) %in% "Data")]

		# check existing input:
		if(!is.data.frame(Data)){
			stop("Argument 'Data' needs to be of class 'data.frame'!")
		}
        if (data.table::is.data.table(Data)) {
            Data <- as.data.frame(Data)
        }
		nms <- names(Data)
		sind_Input <- sNames13 %in% nms
		if(any(sind_Input)){
			sData <- Data[,sNames13[sind_Input],drop=FALSE]
			Data <- Data[,-which(nms %in% sNames13[sind_Input]),drop=FALSE]
		} else {
			sData <- NULL
		}
		nms <- names(Data)
		cind_Input <- cNames13 %in% nms
		if(any(cind_Input)){
			cData <- Data[,cNames13[cind_Input],drop=FALSE]
			Data <- Data[,-which(nms %in% cNames13[cind_Input]),drop=FALSE]
			names(cData) <- sNames13[cind_Input]
		} else {
			cData <- NULL
		}

		# check duplicated names:
		nms <- names(Data)
		if(!is.null(nms)){
			nms0 <- intersect(nms,c(sNames,cNames))
			if(length(nms0)){
				nms1 <- paste0(nms0," (Data)")
				warning("Reserved column names:\n",paste(paste0("-> ",nms0),collapse="\n"),"\nChanging to:\n",paste(paste0("-> ",nms1),collapse="\n"))
				names(Data)[match(nms0,nms)] <- nms1
			}
		}
		# output
		bind <- list(sData,cData,Data)
		Out <- do.call(cbind,bind[!sapply(bind,is.null)])
		nrOut <- NROW(Out)
		if(length(mc) > 0){
			DFmc <- as.data.frame(mc,stringsAsFactors=FALSE)
			nrDF <- NROW(DFmc)			
			if(as.logical(nrOut %% nrDF)){
				stop("Number of Data rows is not a multiple of Argument length(s)!")
			}
			ind <- rep(1:nrDF,nrOut/nrDF)
			nmc <- names(mc)
			snm <- sNames[c(which(sind_Input),which(cind_Input))]
			if(any(snm %in% nmc)){
				nmadd <- intersect(nmc,snm)
				DFrep <- as.data.frame(mc[nmadd],stringsAsFactors=FALSE)
				Out[,nmadd] <- DFrep[ind,]
				nmc <- nmc %w/o% nmadd
			}
			addOut <- as.data.frame(mc[nmc],stringsAsFactors=FALSE)
			names(addOut) <- nmc
			Out <- cbind(Out,addOut[ind,,drop=FALSE])
		}

		# missing args
		margs <- sNames13 %w/o% names(Out)
		if(length(margs)){
			addOut <- as.data.frame(lapply(margs,get,envir=environment()),stringsAsFactors=FALSE)
			names(addOut) <- margs
			Out <- cbind(Out,addOut[rep(1,nrOut),,drop=FALSE])
		}

		# order
		nms_add <- names(Out) %w/o% sNames13
		Out <- Out[,c(sNames13,nms_add)]
		names(Out)[1:13] <- cNames13
		row.names(Out) <- seq_along(Out[,1])
	} else {
		Out <- as.data.frame(lapply(sNames13,get,envir=environment()),stringsAsFactors=FALSE)
		names(Out) <- cNames13
	}

	if(is.factor(Out[,12])) Out[,12] <- as.character(Out[,12])
	if(is.factor(Out[,13])) Out[,13] <- as.character(Out[,13])
	if(any(!sapply(Out[,1:11],is.numeric))|any(!sapply(Out[,c(12,13)],is.character)))stop("input must be of type numeric or character (SensorNames and SourceNames)!\n")
	if(any(Out[,"N0"] < 1E3))warning("N0 should be larger than 1000! (preferably larger than 1E4!)\n")

	Out[,"Ustar [m/s]"] <- round(Out[,"Ustar [m/s]"],7)
	Out[which(Out[,"L [m]"] > 10^9),"L [m]"] <- 10^9
	Out[which(Out[,"L [m]"] < -10^9),"L [m]"] <- -10^9
	Out[,"L [m]"] <- round(Out[,"L [m]"],1)
	Out[,"L [m]"] <- ifelse(abs(Out[,"L [m]"])<0.1,NA,Out[,"L [m]"])
	suppressWarnings(Out[,"Zo [m]"] <- as.numeric(sprintf("%1.3e",Out[,"Zo [m]"])))
	Out[,"SigmaU/Ustar [-]"] <- round(Out[,"SigmaU/Ustar [-]"],3)
	Out[,"SigmaV/Ustar [-]"] <- round(Out[,"SigmaV/Ustar [-]"],3)
	Out[,"SigmaW/Ustar [-]"] <- round(Out[,"SigmaW/Ustar [-]"],3)
	Out[,"SigmaW/Ustar Height [m]"] <- round(Out[,"SigmaW/Ustar Height [m]"],3)
	Out[,"WD [deg N]"] <- round(Out[,"WD [deg N]"],3)
	Out[,"d [m]"] <- round(Out[,"d [m]"],3)
	Out[,"N0"] <- round(Out[,"N0"],0)
	Out[,"MaxFetch"] <- round(Out[,"MaxFetch"],0)

	# check abs(rho_uw) <= 1
	isna <- as.logical(rowSums(is.na(Out[,1:13])))
	if(any(wruw <- abs(1/Out[!isna,"SigmaU/Ustar [-]"]/calcbw(Out[!isna,"SigmaW/Ustar [-]"],(Out[!isna,"SigmaW/Ustar Height [m]"]-Out[!isna,"d [m]"])/Out[!isna,"L [m]"]))>1)){
		cat("Line:\n")
		cat(paste(wruw <- which(wruw),collapse=","),"\n")
		cat("Correlation of u and w is greater than 1! -> Setting Ustar value to NA!\n")
		Out[!isna,][wruw,2] <- NA
	}

	return(structure(Out,class=c("Interval","data.frame")))
}
