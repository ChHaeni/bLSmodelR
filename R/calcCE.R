.calcCE <- function(SubRun, InputList, Srcs, C.Path, variables = 'CE') {

    # which additional variables?
    which_vars <- c('uCE', 'vCE', 'wCE') %in% variables

    # record gc/mem
    .record_now(start = TRUE)

	# SubRun <- SncRun[ilist,]
	# zu beginn on.exit(fehlerangabe!?)
	setDT(SubRun)

    # get number of trajectories
	N0 <- SubRun[1, N0]

    # get all Sensors
	AllSensorNames <- unlist(strsplit(SubRun[, Calc.Sensor], ",", fixed = TRUE))
	lasn <- length(AllSensorNames)

    # get all Sources
	SourceNames <- unlist(strsplit(SubRun[1, Source], ",", fixed = TRUE))	
	Scalc <- copy(Srcs[SourceNames])
	setkey(Scalc, Plot, pid)

	# initialize CE list
	Ci_sub <- vector(mode = "list", length = lasn)
	names(Ci_sub) <- AllSensorNames

    # final CE list per Source
	Ci <- lapply(seq_along(SourceNames), function(x) Ci_sub)
	names(Ci) <- SourceNames

	# initialize uvw list
	UVW <- vector(mode = "list", length = nrow(SubRun))

	# initialize Output
	Out <- SubRun[rep(1, length(SourceNames)), ][, 
        SensorHeight := as.character(SensorHeight)]
    # get Sensor heights
	Sheight <- range(InputList$Sensors$"Calc.Sensors"[
        InputList[['Sensors']][['Calc.Sensors']][, 'Sensor Name'] == Out[, Sensor]
        , "Sensor Height (m)"])
	if (length(Sheight) > 1) {
		Sheight <- paste(sprintf("%1.2f", Sheight), collapse = " to ")
	} else {
		Sheight <- sprintf("%1.3f", Sheight[1])
	}
    # fill & prepare Out
	Out[, ":="(
			Source = SourceNames, 
            SourceArea = attr(InputList[["Sources"]], "SAreas")[SourceNames], 
			SensorHeight = Sheight, 
			CE = 0, CE_se = NA_real_, CE_lo = NA_real_, CE_hi = NA_real_, 
			uCE = 0, uCE_se = NA_real_, uCE_lo = NA_real_, uCE_hi = NA_real_, 
			vCE = 0, vCE_se = NA_real_, vCE_lo = NA_real_, vCE_hi = NA_real_, 
			wCE = 0, wCE_se = NA_real_, wCE_lo = NA_real_, wCE_hi = NA_real_, 
			N_TD = 0, TD_Time_avg = NA_real_, TD_Time_max = NA_real_, 
            Max_Dist = NA_real_, UCE = 0, n_time_avg = 0, Calc.Sensor = NULL, 
            seed = NULL)]
    # remove Calc.* output
    rm_ind <- Reduce(':', match(c('Calc.mtime', 'Calc.N0'), names(Out)))
    Out[, names(Out)[rm_ind] := NULL]
    Out[, Sensor_Swustar := NULL]
	setkey(Out, Source)

	for(Row in seq(nsr <- nrow(SubRun))){

		SensorNames <- unlist(strsplit(SubRun[Row, Calc.Sensor], ",", fixed = TRUE))
		
        # TODO: Srange by row -> same for deposition (by Catalog)
        # prepare range of Source - Sensor distance
        sind <- chmatch(SensorNames, InputList$Sensors$"Calc.Sensors"[, "Point Sensor Name"])
        SensorPositions <- as.matrix(InputList$Sensors$"Calc.Sensors"[sind, c("x-Coord (m)", "y-Coord (m)")])
        rownames(SensorPositions) <- InputList$Sensors$"Calc.Sensors"[sind, "Point Sensor Name"]
        Srange <- Scalc[, rbind(
            cbind(
                x = max(x) - min(SensorPositions[, "x-Coord (m)"]),
                y = max(y) - min(SensorPositions[, "y-Coord (m)"])
            ),
            cbind(
                x = min(x) - max(SensorPositions[, "x-Coord (m)"]),
                y = min(y) - max(SensorPositions[, "y-Coord (m)"])
            )
        )]

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
		tag_bbox(Catalog,Srange)
		Catalog[, inside_Srange := bbox_inside]
		Catalog[,rn:=.I]

		cat("\n~~ Sensor Height (z-d) >",SubRun[Row,SensorHeight],"m < ~~\n")
		if (Catalog[, any(bbox_inside)]) {
			combs <- expand.grid(SourceNames,SensorNames,KEEP.OUT.ATTRS=FALSE,stringsAsFactors=FALSE)
			SensorNumbers <- chmatch(combs[,2],AllSensorNames)
			nc <- NROW(combs)
			Steps <- unique(round(seq(1,nc,length.out=11))) %w/o% nc
			for(cmb in seq(nc)){
				
				# progressbar
				if(cmb %in% Steps)cat(sprintf("\r[%s%s%s] %1.0f%%",paste(rep(">",round(2*(cmb-1)/nc*10)),collapse=""),"|",paste(rep(".",max(0,19-round(2*(cmb-1)/nc*10))),collapse=""),(cmb-1)/nc*100))
				
				Source <- combs[cmb,1]
				Sensor <- combs[cmb,2]
				SourceAreaRelative <- copy(Scalc[Source])[, ":="(
                    x = x - SensorPositions[Sensor, 1], 
                    y = y - SensorPositions[Sensor, 2]
                    )]
				Catalog[, bbox_inside := inside_Srange]
				tag_bbox(Catalog, SourceAreaRelative)
				
				if (Catalog[, any(bbox_inside)]) {
					# tag Inside Source
					TDinside <- SourceAreaRelative[, {
						tag_bbox(Catalog, .(x = x, y = y))
						cbind(
                            ID = Catalog[(bbox_inside), rn], 
                            pnt.in.poly(Catalog[(bbox_inside), cbind(x, y)], cbind(x, y))
                        )
					}, by = pid][, sum(pip), by = ID]

					if (TDinside[, any(V1 > 0)]) {
						setkey(Catalog, rn)
                        Catalog[, td_inside := FALSE]
						Catalog[TDinside, td_inside := V1 > 0]

						# calc CE
						Ci[[Source]][[Sensor]] <- Catalog[(td_inside),
                            .(CE = sum(2 / wTD))
                            , by = Traj_ID]
						# Max_Dist etc. (Max_Dist = max fetch inside Source area)
						setkey(Catalog, Traj_ID)
                        # Time >= minTime excludes TD before first TD inside source
						Cat <- Catalog[Catalog[(td_inside), 
                                .(minTime = min(Time))
                                , by = Traj_ID]][Time >= minTime, ]
                        # sum is used because of na.rm option
                        avg_time <- Catalog[(td_inside), -mean(Time)]
						Out[Source, ":="(
							Max_Dist = max(
                                Max_Dist, 
                                rotateCatalog(Cat, SubRun[Row, WD], 
                                    back = TRUE)[, -min(x)], 
                                na.rm = TRUE), 
							N_TD = N_TD + Catalog[, sum(td_inside)], 
							TD_Time_avg = sum(
                                TD_Time_avg, 
                                avg_time,
                                na.rm = TRUE), 
							TD_Time_max = max(
                                TD_Time_max, 
                                Catalog[(td_inside), -min(Time)], 
                                na.rm = TRUE), 
							n_time_avg = n_time_avg + as.integer(!is.na(avg_time))
							)]
					}
				}
			}
			cat(paste0("\r[",paste0(rep(">",20),collapse=""),"] 100%\n"))
		}
        # record gc/mem
        .record_now()
	}
	rm(Catalog)

	cat("\nCalculating Source contributions\n\n")

	# TD_Time_avg
	if (Out[, any(N_TD > 0)]) {
        Out[N_TD > 0, TD_Time_avg := TD_Time_avg / n_time_avg]
    }
	Out[, n_time_avg := NULL]

	# weights:
	wts <- rep(2,lasn)
	wts[c(1,lasn)] <- 1
	rwts <- wts/sum(wts)
	if(lasn > 1){
		# rwts korrekt sortieren -> welche Sensoren sind max/min node? -> welche pos haben diese?
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
    UVW_mean <- lapply(UVW,function(x)colMeans(x))

    if (which_vars[1]) {
        # U Matrix:
        U_matrix <- matrix(sapply(1:nrow(SubRun),function(x)UVW[[x]][,"u0"] - UVW_mean[[x]]["u0"]),nrow=N0)
    }
    if (which_vars[2]) {
        V_matrix <- matrix(sapply(1:nrow(SubRun),function(x)UVW[[x]][,"v0"] - UVW_mean[[x]]["v0"]),nrow=N0)
    }
    if (which_vars[3]) {
        W_matrix <- matrix(sapply(1:nrow(SubRun),function(x)UVW[[x]][,"w0"] - UVW_mean[[x]]["w0"]),nrow=N0)
    }

	for(i in Out[N_TD>0,Source]){

		# check NULL
		not_null <- !sapply(Ci[[i]], is.null)
		
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
        if (which_vars[1]) {
            uvwCE <- c_matrix*U_matrix[,uvw_key[AllSensorNames,row]]
            uCE_add <- sum(colSums(uvwCE)*rwts)/N0
            uCE_se_add <- sqrt(sum(cov(uvwCE)*orwts)/N0)
        } else {
            uCE_add <- uCE_se_add <- NA_real_
        }

		# vCE + SE
        if (which_vars[2]) {
            uvwCE <- c_matrix*V_matrix[,uvw_key[AllSensorNames,row]]
            vCE_add <- sum(colSums(uvwCE)*rwts)/N0
            vCE_se_add <- sqrt(sum(cov(uvwCE)*orwts)/N0)
        } else {
            vCE_add <- vCE_se_add <- NA_real_
        }

		# uCE + SE
        if (which_vars[3]) {
            uvwCE <- c_matrix*W_matrix[,uvw_key[AllSensorNames,row]]
            wCE_add <- sum(colSums(uvwCE)*rwts)/N0
            wCE_se_add <- sqrt(sum(cov(uvwCE)*orwts)/N0)
        } else {
            wCE_add <- wCE_se_add <- NA_real_
        }

		# write results
		Out[i, CE := CE_mean]
		Out[i, UCE := UCE_mean]
		Out[i, uCE := uCE_add]
		Out[i, vCE := vCE_add]
		Out[i, wCE := wCE_add]
		# SE:
		Out[i, CE_se := CE_se_add]
		Out[i, uCE_se := uCE_se_add]
		Out[i, vCE_se := vCE_se_add]
		Out[i, wCE_se := wCE_se_add]

	}
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
    # record gc/mem
    cpu_mem <- .record_now(reset = TRUE)
    setattr(Out, 'cpu_mem', cpu_mem)

	return(Out)

}
