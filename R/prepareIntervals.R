prepareIntervals <- function(InputList,C.Path=NULL,asDT=TRUE,simpleNames=TRUE,ncores=1){
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
	Tol[,which(Tol[,1:6] == 0)] <- 1E-3


	# create extended Cat_Name name
	IntExt[,Cat.Name := {
		gsub("_$",sprintf("-%05.0f_%05.0f_%05.0f",.BY$sUu*1E4,.BY$sVu*1E4,.BY$Sensor_Swustar*1E4),createCatName(.SD))
	},by=.(sUu,sVu,Sensor_Swustar)]


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

		if(nrow(CList) > 0){
			# rename
			setnames(CList,names(CList)[-1],paste0("Cat_",names(CList)[-1]))
			# create Cat_Sensor_Swustar
			CList[,Cat_Sensor_Swustar := round(calcsigmaW(1,Cat_ZSens/Cat_L,Cat_bw),3)]

			if(!is.null(cl) && nrow(IntExt) >= length(cl)){
				# t1p <- Sys.time()
				# split index
                ind <- parallel::clusterSplit(cl, seq.int(nrow(IntExt)))
				Key <- rbindlist(parallel::clusterApply(cl, lapply(ind, function(x, y)y[x, ], y = IntExt[, c(.SD, .(
					sUu_Upper = sUu*TolUpper[,"SigmaU/Ustar"],
					sUu_Lower = sUu*TolLower[,"SigmaU/Ustar"],
					sVu_Upper = sVu*TolUpper[,"SigmaV/Ustar"],
					sVu_Lower = sVu*TolLower[,"SigmaV/Ustar"],
					sWu_Upper = Sensor_Swustar*TolUpper[,"SigmaW/Ustar"],
					sWu_Lower = Sensor_Swustar*TolLower[,"SigmaW/Ustar"],
					Zo_Upper = Zo*TolUpper[,"Zo"],
					Zo_Lower = Zo*TolLower[,"Zo"],
					L_Upper = abs(L)*TolUpper[,"L"],
					L_Lower = abs(L)*TolLower[,"L"],
					SensorHeight_Upper = SensorHeight*TolUpper[,"Sensor Height"],
					SensorHeight_Lower = SensorHeight*TolLower[,"Sensor Height"]
					))]),.MatchWrapper,Cat_list=CList,Tol_=Tol))
				# t2p <- Sys.time()
			} else { # not parallel

				# t1s <- Sys.time()

				# check matching
				Key <- IntExt[,c(.SD,.(
					sUu_Upper = sUu*TolUpper[,"SigmaU/Ustar"],
					sUu_Lower = sUu*TolLower[,"SigmaU/Ustar"],
					sVu_Upper = sVu*TolUpper[,"SigmaV/Ustar"],
					sVu_Lower = sVu*TolLower[,"SigmaV/Ustar"],
					sWu_Upper = Sensor_Swustar*TolUpper[,"SigmaW/Ustar"],
					sWu_Lower = Sensor_Swustar*TolLower[,"SigmaW/Ustar"],
					Zo_Upper = Zo*TolUpper[,"Zo"],
					Zo_Lower = Zo*TolLower[,"Zo"],
					L_Upper = abs(L)*TolUpper[,"L"],
					L_Lower = abs(L)*TolLower[,"L"],
					SensorHeight_Upper = SensorHeight*TolUpper[,"Sensor Height"],
					SensorHeight_Lower = SensorHeight*TolLower[,"Sensor Height"]
					))][,{
					CList[
						abs(Cat_kv - .BY$kv) < 1E-2 &
						abs(Cat_A - .BY$A) < 1E-2 &
						abs(Cat_alpha - .BY$alpha) < 1E-3 &
						Cat_MaxFetch >= .BY$MaxFetch &
						Cat_ZSens <= SensorHeight_Upper &
						Cat_ZSens >= SensorHeight_Lower &
						sign(Cat_L) == sign(.BY$L) &
						abs(Cat_L) <= L_Upper &
						abs(Cat_L) >= L_Lower &
						Cat_Zo <= Zo_Upper &
						Cat_Zo >= Zo_Lower &
						Cat_Su_Ustar <= sUu_Upper &
						Cat_Su_Ustar >= sUu_Lower &
						Cat_Sv_Ustar <= sVu_Upper &
						Cat_Sv_Ustar >= sVu_Lower &
						Cat_Sensor_Swustar <= sWu_Upper &
						Cat_Sensor_Swustar >= sWu_Lower
					,.SD]
					},by = .(rn,Cat.Name,SensorHeight,L,Zo,sUu,sVu,Sensor_Swustar,MaxFetch,Sensor,N0,alpha,A,kv)][
						,":="(
						devZSens = abs(Cat_ZSens/SensorHeight - 1)/Tol[1],
						devL = abs(Cat_L/L - 1)/Tol[2],
						devZo = abs(Cat_Zo/Zo - 1)/Tol[3],
						devsUu = abs(Cat_Su_Ustar/sUu - 1)/Tol[4],
						devsVu = abs(Cat_Sv_Ustar/sVu - 1)/Tol[5],
						devSensor_Swustar = abs(Cat_Sensor_Swustar/Sensor_Swustar - 1)/Tol[6],
						devN0 = Cat_N0 - N0
						)
					][
						devZSens <= 1.0000001 &
						devL <= 1.0000001 &
						devZo <= 1.0000001 &
						devsUu <= 1.0000001 &
						devsVu <= 1.0000001 &
						devSensor_Swustar <= 1.0000001,sumDev:=devZSens+devL+devZo+devsUu+devsVu+devSensor_Swustar]

				# t2s <- Sys.time()

			} # end if(parl)

			if(nrow(Key) > 0){
				# note to myself: apply biases in sumDev someday
				Key[devN0>0,sumDev:=sumDev+0.001]
				Key[devN0<0,sumDev:=sumDev+6+log(-devN0)]
				
				# get best matches
				Key <- Key[,.(Cat_Name=Name[which.min(sumDev)],Cat_N0=Cat_N0[which.min(sumDev)],R_N0=N0[which.min(sumDev)]
					,R_Name=gsub("-[0-9]{5,}_[0-9]{5,}_[0-9]{5,}$","",Cat.Name[which.min(sumDev)])),by=.(R_Sensor=Sensor,R_rn=rn)]

				# check max N0 per Catalog & define Sensor/rn to calculate
				setkey(CList,Name)
				Key[,c(names(CList)[-1],"K_exists","K_calc","K_N0") := {
					cbind(
						CList[.BY,-1]
						,TRUE
						,c(CList[.BY,Cat_N0] < max(R_N0),rep(FALSE,.N-1))
						,max(R_N0)
					)},by = Cat_Name]
				setkey(Key,R_Sensor,R_rn)	

				rm(CList)

				# merge with IntExt
				IntExt[Key,c(gsub("^Cat_","Calc.",grep("^Cat_",names(Key),value=TRUE)[-(1:2)]),"Cat.Name","Cat.exists","Cat.calc","Calc.N0") := {
					c(mget(grep("^Cat_",names(Key),value=TRUE)[-(1:2)]),list(Cat_Name,K_exists,K_calc,K_N0))
				}]

				cat("** -> Found",Key[,length(unique(Cat_Name))],"matching catalogs.",if(Key[,sum(K_calc)]) paste0("(",Key[,sum(K_calc)]," need extending)\n") else "\n")
				rm(Key)
			} else {
				rm(CList)
				cat("** -> No matching catalog found.\n")
			}
		} else {
			cat("** No TD catalogs to check.\n")
		}


		if(IntExt[,any(!(Cat.exists))]){
			cat("* Checking for cross-matches in supplied intervals...\n")
			
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
				Cat_Sensor = Sensor
				)])
			CatList[,Zeile := seq.int(.N)]
			setkey(CatList,Name)

			if(!is.null(cl) && IntExt[,sum(!Cat.exists)] >= length(cl)){
				# t1p <- Sys.time()
				# split index
                ind <- parallel::clusterSplit(cl, seq.int(IntExt[,sum(!Cat.exists)]))
				Key <- rbindlist(parallel::clusterApply(cl,
                        lapply(ind, function(x, y)y[x, ], y = IntExt[!(Cat.exists), 
                            .(rn, Cat.Name, SensorHeight, L, Zo, sUu, sVu, Sensor_Swustar, MaxFetch, Sensor, N0, 
                                sUu_Upper = sUu*TolUpper[,"SigmaU/Ustar"],
                                sUu_Lower = sUu*TolLower[,"SigmaU/Ustar"],
                                sVu_Upper = sVu*TolUpper[,"SigmaV/Ustar"],
                                sVu_Lower = sVu*TolLower[,"SigmaV/Ustar"],
                                sWu_Upper = Sensor_Swustar*TolUpper[,"SigmaW/Ustar"],
                                sWu_Lower = Sensor_Swustar*TolLower[,"SigmaW/Ustar"],
                                Zo_Upper = Zo*TolUpper[,"Zo"],
                                Zo_Lower = Zo*TolLower[,"Zo"],
                                L_Upper = abs(L)*TolUpper[,"L"],
                                L_Lower = abs(L)*TolLower[,"L"],
                                SensorHeight_Upper = SensorHeight*TolUpper[,"Sensor Height"],
                                SensorHeight_Lower = SensorHeight*TolLower[,"Sensor Height"]
					)]),.CrossMatchWrapper,Cat_list=CatList,Tol_=Tol))
				# t2p <- Sys.time()
			} else {
	      # not parallel
	      
	      # t1s <- Sys.time()

	      # check matching heights etc.
	      CheckPart1 <- IntExt[!(Cat.exists), .(rn, SensorHeight, Sensor_Swustar, MaxFetch,
	        sWu_Upper = Sensor_Swustar*TolUpper[,"SigmaW/Ustar"],
	        sWu_Lower = Sensor_Swustar*TolLower[,"SigmaW/Ustar"],
	        SensorHeight_Upper = SensorHeight*TolUpper[,"Sensor Height"],
	        SensorHeight_Lower = SensorHeight*TolLower[,"Sensor Height"]
	        )][,{
	        CatList[
	          Cat_MaxFetch >= MaxFetch[1] &
	          Cat_Sensor_Swustar <= sWu_Upper[1] &
	          Cat_Sensor_Swustar >= sWu_Lower[1] &
	          Cat_ZSens <= SensorHeight_Upper[1] &
	          Cat_ZSens >= SensorHeight_Lower[1],{
	          .(
	            SensorHeight = SensorHeight[1], 
	            Sensor_Swustar = Sensor_Swustar[1], 
	            MaxFetch = MaxFetch[1], 
	            Zeile = Zeile)
	        }]    
	      }, keyby = .(ssm = paste(SensorHeight, Sensor_Swustar, MaxFetch, sep = "/"))]

	      # check matching MOST
	      Key <- IntExt[!(Cat.exists),.(rn,Cat.Name,SensorHeight,L,Zo,sUu,sVu,Sensor_Swustar,MaxFetch,Sensor,N0,
	        sUu_Upper = sUu*TolUpper[,"SigmaU/Ustar"],
	        sUu_Lower = sUu*TolLower[,"SigmaU/Ustar"],
	        sVu_Upper = sVu*TolUpper[,"SigmaV/Ustar"],
	        sVu_Lower = sVu*TolLower[,"SigmaV/Ustar"],
	        Zo_Upper = Zo*TolUpper[,"Zo"],
	        Zo_Lower = Zo*TolLower[,"Zo"],
	        L_Upper = abs(L)*TolUpper[,"L"],
	        L_Lower = abs(L)*TolLower[,"L"]
	        )][,{

	        ind <- paste(SensorHeight, Sensor_Swustar, MaxFetch, sep = "/")
	        Sub <- CheckPart1[.(ind), .(Zeile, ssm)]
	        ind2 <- Sub[,unique(Zeile)]

	        out <- CatList[match(ind2,Zeile)][
	          Cat_Su_Ustar <= sUu_Upper[1] &
	          Cat_Su_Ustar >= sUu_Lower[1] &
	          Cat_Sv_Ustar <= sVu_Upper[1] &
	          Cat_Sv_Ustar >= sVu_Lower[1] &
	          Cat_Zo <= Zo_Upper[1] &
	          Cat_Zo >= Zo_Lower[1] &
	          sign(Cat_L) == sign(L[1]) &
	          abs(Cat_L) <= L_Upper[1] &
	          abs(Cat_L) >= L_Lower[1],
	            merge(.SD, Sub, by = "Zeile")]

	        c(
	          .SD[out[,match(ssm, ind)],],
	          out
	          )
	      }, by = rn][,":="(
		        devZSens = abs(Cat_ZSens/SensorHeight - 1)/Tol[1],
		        devL = abs(Cat_L/L - 1)/Tol[2],
		        devZo = abs(Cat_Zo/Zo - 1)/Tol[3],
		        devsUu = abs(Cat_Su_Ustar/sUu - 1)/Tol[4],
		        devsVu = abs(Cat_Sv_Ustar/sVu - 1)/Tol[5],
		        devSensor_Swustar = abs(Cat_Sensor_Swustar/Sensor_Swustar - 1)/Tol[6],
		        devN0 = Cat_N0 - N0)][
	        devZSens <= 1.0000001 &
	        devL <= 1.0000001 &
	        devZo <= 1.0000001 &
	        devsUu <= 1.0000001 &
	        devsVu <= 1.0000001 &
	        devSensor_Swustar <= 1.0000001,sumDev:=devZSens+devL+devZo+devsUu+devsVu+devSensor_Swustar
	        ]

	      # t2s <- Sys.time()
	      rm(CheckPart1)
	      #   t2s - t1s

	        

			}

			if(nrow(Key)==IntExt[,sum(!Cat.exists)]){
				cat("** No cross-matching rows...\n")
				IntExt[!(Cat.exists),Cat.calc := TRUE]
				# assign Calc.values
				IntExt[,":="(
					Calc.ZSens=SensorHeight,
					Calc.L=L,
					Calc.Zo=Zo,
					Calc.Su_Ustar=sUu,
					Calc.Sv_Ustar=sVu,
					Calc.Sensor_Swustar=round(calcsigmaW(1,SensorHeight/L,bw),3),
					Calc.bw=bw,
					Calc.MaxFetch=MaxFetch,
					Calc.C0 = C0,
					Calc.N0 = N0,
					Calc.Ustar = Ustar,
					Calc.alpha = alpha,
					Calc.A = A,
					Calc.kv = kv
					)]
			} else {
				# note to myself: apply biases in sumDev someday
				Key[devN0>0,sumDev:=sumDev+0.001]
				Key[devN0<0,sumDev:=sumDev+6+log(-devN0)]

				# check for most (& best) matching
				cat("** Checking best cross-matching catalogs...\n")
				setorder(Key[,N:=.N,by=Name],-N, sumDev)
				keys <- Key[,unique(Name)]
				setkey(Key,Sensor,rn)

				CatNames <- gsub("^Cat_","Calc.",grep("^Cat_",names(CatList)[2:10],value=TRUE))
				nk <- length(keys)
				nk10 <- max(round(nk/10), 1)
				cat(paste0(sprintf("\r** %3.0f%% [",0),paste(rep(c(">>","  "),c(0,10 - 0)),collapse=""),"]"))
				for(i in seq_along(keys)){
					if(i%%nk10 == 0){
						np <- i*100/nk
						np10 <- round(np/10)
						cat(paste0(sprintf("\r** %3.0f%% [",np),paste(rep(c(">>","  "),c(np10,10 - np10)),collapse=""),"]"))
					}
					k <- keys[i]
					if(k != ""){
						k1 <- Key[Name %chin% k, .(Sensor, rn, Cat.Name)]
						IntExt[k1,c("Cat.Name","Calc.N0","Calc.Ustar","Calc.alpha","Calc.A","Calc.kv",CatNames) := {
							c(list(k,max(N0),Ustar,alpha,A,kv),CatList[.(k),names(CatList)[2:10],with=FALSE][which.max(Cat_MaxFetch)])
						}]
						# set inheriting rows with different Cat.Name
						IntExt[k1[!(Cat.Name %chin% k)],Cat.calc := FALSE]
						# set dominant row + inheriting rows with identical Cat.Name
						IntExt[(Cat.Name %chin% k),Cat.calc := {
							out <- rep(FALSE,.N)
							out[which.max(MaxFetch)] <- TRUE
							out
						}]
						### 'delete' rows in Key
						# keys
						keys[keys %chin% k1[,Cat.Name]] <- ""
						# k1 rn/Sensor
						Key[k1,Name := ""]
						# rows depending on k1 rows
						Key[Name %chin% k1[,Cat.Name],Name := ""]
					}
				}
				cat("\n** done.\n")
			}

			rm(Key,CatList)
		}
	} else {
		# TDread == FALSE: 
		IntExt[,":="(
			Calc.ZSens=SensorHeight,
			Calc.L=L,
			Calc.Zo=Zo,
			Calc.Su_Ustar=sUu,
			Calc.Sv_Ustar=sVu,
			Calc.Sensor_Swustar=round(calcsigmaW(1,SensorHeight/L,bw),3),
			Calc.bw=bw,
			Calc.MaxFetch=MaxFetch,
			Calc.C0 = C0,
			Calc.N0 = N0,
			Calc.Ustar = Ustar,
			Calc.alpha = alpha,
			Calc.A = A,
			Calc.kv = kv
			)]
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
