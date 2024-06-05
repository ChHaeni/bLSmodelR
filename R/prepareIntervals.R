prepareIntervals <- function(InputList, C.Path = NULL, asDT = TRUE, simpleNames = TRUE, 
    ncores = 1, throttle = 100, skipCrossCheck = FALSE) {
    # throttle argument idea taken from ?setDTthreads. Value is arbitrary
	cat("Preparing intervals for model run:\n")
	# check NA in Interval data.frame:
	whichNA <- as.logical(rowSums(is.na(InputList[["Interval"]][,1:13])))
	if(any(whichNA)){
		cat(paste0("* Cleaning Interval data.frame -> Removed ",sum(whichNA)," row",if(sum(whichNA)>1) "s"," containing NA values.\n"))
	}
	
	#### extend Interval:
	cat("* Restructuring intervals...\n")
	# create data.table:
	Int <- data.table(InputList[["Interval"]][!whichNA,],keep.rownames=TRUE)
	# set row order
	Int[, orignal__order := seq_len(.N)]

	if((lsnc <- length(Int))>14)setnames(Int,names(Int)[15:lsnc],paste0(names(Int)[15:lsnc],"_removeMeLater"))
	setnames(Int,oldNames <- names(InputList[["Interval"]])[1:13],c("Ustar","L","Zo","sUu","sVu","sWu","z_sWu","WD","d","N0","MaxFetch","Sensor","Source"))

	### procSensors:
	pSens <- procSensors(InputList[["Sensors"]])

	# update Sensor Names (negative lookahead Buchstabe, positive lookbehind Komma oder Anfang)
	for(sn in names(pSens$PS_list)){
		Int[, Sensor := gsub(
			paste0("\\b", sn, "\\b"),
			pSens$PS_list[[sn]], 
			Sensor)]
	}

	# add additional info
	Int[,":="(Cat.Name="",Cat.exists=FALSE,Cat.calc=TRUE,kv=InputList[["Model"]][["kv"]],A=InputList[["Model"]][["A"]],alpha=InputList[["Model"]][["alpha"]])]
	
	Int[,":="(z_sWu = z_sWu - d,
				SensorHeight = 0.0)][
			,bw := round(calcbw(sWu,z_sWu/L),3)
			][
				,":="(C0 = round(calcC0(bw,InputList[["Model"]][["kv"]],InputList[["Model"]][["A"]]),3),
					Sensor_Swustar = 0.0)
			]

	# individual heights
	nhts <- pSens$"heights"
	IntExt <- Int[,{
		# get all Sensors
		iSens <- unlist(strsplit(.BY$Sensor,split=","))
		# index by height
		posL <- sapply(iSens,
			function(sn){
				which(sapply(nhts, function(x, y) any(x %in% y), y = sn, USE.NAMES = FALSE))
			}, USE.NAMES = FALSE)
		# paste by height
		Snames <- tapply(iSens, posL, paste, collapse = ",")
		# write out
		n <- .N
		.SD[rep(seq.int(n),each=length(Snames)),][,":="(
			Newname = rep(Snames,n)
			,SensorHeight = rep(as.numeric(names(nhts))[as.numeric(names(Snames))], n)
			)]
	},by = Sensor][,Sensor := NULL]
	setnames(IntExt,"Newname","Sensor")
	setcolorder(IntExt,names(Int))
	rm(Int)
	
	IntExt[,Sensor_Swustar := round(calcsigmaW(1,SensorHeight/L,bw),3)]
	setkey(IntExt,Sensor,rn)
	####

    cl <- NULL
    if (is.null(ncores)) {
        ncores <- 1
    } else if (inherits(ncores, 'cluster')) {
        cl <- ncores
        ncores <- length(cl)
    } else if (ncores > 1) {
		on.exit(parallel::stopCluster(cl))
        cl <- parallel::makePSOCKcluster(ncores)
    } else if (ncores != 1) {
		stop("Number of cores must be greater or equal to 1!")
    }

	# optimize MaxFetch
	if(IntExt[,any(MaxFetch < 0)]){
		cat("* Optimizing 'MaxFetch'...\n")
		if(is.null(InputList[["Sources"]]))stop("No sources supplied, can not optimize MaxFetch!")

		if(!is.null(cl) && IntExt[,sum(MaxFetch < 0) > ncores]){
            ind <- parallel::clusterSplit(cl, IntExt[, which(MaxFetch < 0)])
			mf <- rbindlist(parallel::clusterApply(cl, lapply(ind, function(x, y)y[x, ], y = IntExt), 
                    .MaxFetchWrapper, p_Sens = pSens, Input_List = InputList))[, MaxFetch]
			IntExt[MaxFetch < 0, MaxFetch := mf]
		} else {
			IntExt[MaxFetch < 0,
				MaxFetch := {
					nmSens <- unlist(strsplit(Sensor,split=","))
					nmSous <- unlist(strsplit(Source,split=","))
					Sens <- structure(pSens$Calc.Sensors[pSens$Calc.Sensors[, "Point Sensor Name"] %in% nmSens, 
						c("x-Coord (m)", "y-Coord (m)")],class="data.frame")
					Sous <- structure(InputList[["Sources"]][InputList[["Sources"]][,1] %in% nmSous,2:3],class="data.frame")
					out <- numeric(length(WD))
					for(i in seq_along(WD)){
						Sensrot <- rotate(Sens,-WD[i])
						Sourot <- rotate(Sous,-WD[i])
						out[i] <- max(Sensrot[,1]) - min(Sourot[,1])  
					}
					out <- ceiling(pmax(out,0)) - MaxFetch
					out
				}
			,by=.(Sensor,Source)]	
		}
	}

	Tol <- t(InputList[["Tolerances"]]/100)

	TolUpper <- 1 + Tol
	TolLower <- 1 - Tol

	# set to very low (because of dividing)
	Tol[,which(Tol[,1:6] == 0)] <- 1E-6


	# create extended Cat_Name name
	IntExt[,Cat.Name := {
		gsub("_$",sprintf("-%05.0f_%05.0f_%05.0f",.BY$sUu*1E4,.BY$sVu*1E4,.BY$Sensor_Swustar*1E4),createCatName(.SD))
	},by=.(sUu,sVu,Sensor_Swustar)]

    # prepare Calc.* columns
    IntExt[, ":="(
        Calc.ZSens = SensorHeight, 
        Calc.L = L, 
        Calc.Zo = Zo, 
        Calc.Su_Ustar = sUu, 
        Calc.Sv_Ustar = sVu, 
        Calc.Sensor_Swustar = round(calcsigmaW(1, SensorHeight / L, bw), 3), 
        Calc.bw = bw, 
        Calc.MaxFetch = MaxFetch, 
        Calc.C0 = C0, 
        Calc.N0 = N0, 
        Calc.Ustar = Ustar, 
        Calc.alpha = alpha, 
        Calc.A = A, 
        Calc.kv = kv
        )]

	#### check existing Catalogs:
	if(InputList[["Model"]][["TDread"]]){

		cat("* Checking for matches in existing TD catalogs...\n")

		### get ranges for catalog testing
		Tol_nms <- paste0("Tol_",c("SensorHeight","L","Zo","sUu","sVu","MaxFetch","alpha","A","kv"))
		Tol_check <- IntExt[,lapply(.(SensorHeight,abs(L),Zo,sUu,sVu,MaxFetch,alpha,A,kv),range)]
		setnames(Tol_check,Tol_nms) 
		Tol_check[1,c("Tol_SensorHeight","Tol_L","Tol_Zo","Tol_sUu","Tol_sVu") := {
			as.list(c(Tol_SensorHeight,Tol_L,Tol_Zo,Tol_sUu,Tol_sVu)*TolLower[1:5])
		}]
		Tol_check[2,c("Tol_SensorHeight","Tol_L","Tol_Zo","Tol_sUu","Tol_sVu") := {
			as.list(c(Tol_SensorHeight,Tol_L,Tol_Zo,Tol_sUu,Tol_sVu)*TolUpper[1:5])
		}]


		if(is.null(C.Path))stop("Please supply a path to a TD catalog folder (argument 'C.Path')!")
		CatList <- data.table(rebuildCatListFile(C.Path))

		# pre-check CatList
		CList <- CatList[
			MaxFetch >= Tol_check[ ,Tol_MaxFetch[1]] &
			alpha == Tol_check[ ,Tol_alpha[1]] &
			kv == Tol_check[ ,Tol_kv[1]] &
			A == Tol_check[ ,Tol_A[1]] & 
			ZSens >= Tol_check[ ,Tol_SensorHeight[1]] & ZSens <= Tol_check[ ,Tol_SensorHeight[2]] &
			abs(L) >= Tol_check[ ,Tol_L[1]] & abs(L) <= Tol_check[ ,Tol_L[2]] &
			Zo >= Tol_check[ ,Tol_Zo[1]] & Zo <= Tol_check[ ,Tol_Zo[2]] &
			Su_Ustar >= Tol_check[ ,Tol_sUu[1]] & Su_Ustar <= Tol_check[ ,Tol_sUu[2]] &
			Sv_Ustar >= Tol_check[ ,Tol_sVu[1]] & Sv_Ustar <= Tol_check[ ,Tol_sVu[2]]
		]
		rm(CatList)

        # add row index
        IntExt[, row := .I]

		if (nrow(CList) > 0) {
            # throttle idea taken from ?setDTthreads
            if (ncores > 1 && IntExt[, uniqueN(Cat.Name) > throttle * ncores]) {
                # split by cat name
                IntSplit <- split(IntExt, by = 'Cat.Name')
                # run load balanced
                out <- .clusterApplyLB(cl, IntSplit, .CheckCatMatches, cat_list = CList,
                    tol = Tol, tol_lower = TolLower, tol_upper = TolUpper)
                # rbind results
                IntExt <- rbindlist(out, fill = TRUE)
                # fix order
                setorder(IntExt, row)
            } else {
                .CheckCatMatches(IntExt, CList, Tol, TolLower, TolUpper)
            }
            # be verbose
            if ('Cat.extend' %in% names(IntExt)) {
                # fix any NA due to parallelism
                IntExt[is.na(Cat.extend), Cat.extend := FALSE]
                # be verbose
                cat("** -> Found", IntExt[, sum(Cat.exists)], "rows matching", 
                    IntExt[(Cat.exists), uniqueN(Cat.Name)], "catalogs.",
                    if (IntExt[, any(Cat.extend)]) paste0("(", IntExt[(Cat.extend), uniqueN(Cat.Name)], 
                        " need extending)\n") else "\n")
                # remove again
                IntExt[, Cat.extend := NULL]
            } else {
                cat("** -> No matching catalog found.\n")
            }
		} else {
			cat("** No TD catalogs to check.\n")
		}

		if (IntExt[, any(!(Cat.exists))]) {
			
			#### check Catalog joins:
			CatList <- unique(IntExt[!(Cat.exists),.(
				Name=Cat.Name,
				Cat_ZSens=SensorHeight,
				Cat_L=L,
				Cat_Zo=Zo,
				Cat_Su_Ustar=sUu,
				Cat_Sv_Ustar=sVu,
				Cat_Sensor_Swustar=round(calcsigmaW(1,SensorHeight/L,bw),3),
				Cat_bw=bw,
				Cat_MaxFetch=MaxFetch,
				Cat_C0 = C0,
				Cat_N0 = N0,
				Cat_rn = rn,
				Cat_Sensor = Sensor,
                Zeile = row
				)])
			setkey(CatList,Name)

            if (skipCrossCheck || sum(Tol[1, ]) <= 1e-5) {
                cat("* Get exact matches in supplied intervals...\n")
                # zero tolerance
                IntExt[!(Cat.exists), Cat.calc := {
                    out <- rep(FALSE, .N)
                    # check MaxFetch
                    out[which.max(MaxFetch)] <- TRUE
                    out
                }, by = .(SensorHeight, Sensor_Swustar, L, Zo, sUu, sVu, kv, A, alpha)]
			} else {
                cat("* Checking for cross-matches in supplied intervals...\n")
                # get calc before for number printing
                n_before <- IntExt[, sum(Cat.calc)]
                if (ncores > 1 && IntExt[!(Cat.exists), uniqueN(Cat.Name) > throttle * ncores]) {
                    # split by cat name
                    IntSplit <- split(IntExt[!(Cat.exists)], by = 'Cat.Name')
                    # run load balanced
                    out <- .clusterApplyLB(cl, IntSplit, .CheckCrossMatches, cat_list = CatList,
                        tol = Tol, tol_lower = TolLower, tol_upper = TolUpper)
                    # rbind results
                    IntExt <- rbindlist(out, fill = TRUE)
                    # fix order
                    setorder(IntExt, row)
                } else {
                    .CheckCrossMatches(IntExt, CatList, Tol, TolLower, TolUpper)
                }
                # be verbose
                if ('cat_calc' %in% names(IntExt)) {
                    # remove helper column
                    IntExt[, cat_calc := NULL]
                    cat("** Done. --> Found", n_before - IntExt[!(Cat.exists), sum(Cat.calc)], "cross-matches",
                        paste0("(",IntExt[!Cat.exists & Cat.calc, uniqueN(Cat.Name)], " catalogs)."),"\n\n")
                } else {
                    cat("** No cross-matching rows...\n")
                }
			}

            rm(CatList)

		}

        # remove helper index column
        IntExt[, row := NULL]

        # assign Calc.values
        IntExt[!(Cat.exists), ':='(
            Calc.ZSens = SensorHeight[Cat.calc], 
            Calc.L = L[Cat.calc], 
            Calc.Zo = Zo[Cat.calc], 
            Calc.Su_Ustar = sUu[Cat.calc], 
            Calc.Sv_Ustar = sVu[Cat.calc], 
            Calc.Sensor_Swustar = round(calcsigmaW(1, SensorHeight[Cat.calc] / L[Cat.calc], bw[Cat.calc]), 3), 
            Calc.bw = bw[Cat.calc], 
            Calc.MaxFetch = MaxFetch[Cat.calc], 
            Calc.C0 = C0[Cat.calc], 
            Calc.N0 = N0[Cat.calc], 
            Calc.Ustar = Ustar[Cat.calc], 
            Calc.alpha = alpha[Cat.calc], 
            Calc.A = A[Cat.calc], 
            Calc.kv = kv[Cat.calc]
        ), by = Cat.Name]

	}	
	cat("Finishing preparation...\n")

	# Melt rn/Cat.Name and order by rn:
	cnames <- names(IntExt)

	# collapse rows
	IntExt <- setcolorder(IntExt[,{
		SNames <- paste0(Sensor,collapse=",")
		.SD[max(1,which(Cat.calc)),][,c("Sensor","SensorHeight","Cat.calc") := .(SNames,Calc.ZSens,Cat.calc)]
	},by=.(rn,Cat.Name)],cnames)

	# fix row order
	IntExt <- IntExt[order(orignal__order_removeMeLater)][, orignal__order_removeMeLater := NULL]

	if(!simpleNames){
		setnames(IntExt,2:14,oldNames)
		setnames(IntExt,c("rn","SigmaW/Ustar Height [m]","SensorHeight","Sensor_Swustar"),c("Original Sonic Row","SigmaW/Ustar Height [m above d]","Sensor Height [m above d]","SigmaW/Ustar @ Sensor Height"))
	}

	if(any(grepl("_removeMeLater$",names(IntExt))))setnames(IntExt,gsub("_removeMeLater$","",names(IntExt)))

	if(!asDT){
		setDF(IntExt)
	}
	cat("Done. --> Need to calculate",IntExt[,sum(Cat.calc)],"new catalogs in total.\n")

	return(IntExt)
}
