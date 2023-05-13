.calcDep_Wrapper <- function(RunElement, spatial = FALSE) {
	setDT(RunElement)
	setkey(RunElement, rn, Sensor)
    out <- .calcDep(
        RunElement,
        get('Catalogs', envir = .GlobalEnv), 
        get('Cat.Path', envir = .GlobalEnv), 
        get('Sources', envir = .GlobalEnv), 
        get('Sensors', envir = .GlobalEnv), 
        get('vDep', envir = .GlobalEnv), 
        get('vDepSpatial', envir = .GlobalEnv),
        is_spatial = spatial,
        ci_fun = if (spatial) fill_Ci_spatial else fill_Ci_homogeneous
    )
    # update vd index
	out[, vd_index := RunElement[, vd_index]]
    out
}

.calcDep_Wrapper_noexport <- function(RunElement, 
        Catalogs,
        Cat.Path,
        Sources,
        Sensors,
        vDep,
        vDepSpatial,
    spatial = FALSE) {
	setDT(RunElement)
	setkey(RunElement, rn, Sensor)
    out <- .calcDep(
        RunElement,
        Catalogs,
        Cat.Path,
        Sources,
        Sensors,
        vDep,
        vDepSpatial,
        is_spatial = spatial,
        ci_fun = if (spatial) fill_Ci_spatial else fill_Ci_homogeneous
    )
    # update vd index
	out[, vd_index := RunElement[, vd_index]]
    out
}

.calcDep <- function(Run, Catalogs, C.Path, Sources, CSnsrs, vd, vdSpatial,
    is_spatial, ci_fun) {

    # record gc/mem
    .record_now(start = TRUE)

	# on.exit(browser())
	vd_index <- Run[, vd_index]
	# browser()	
	if(Run[,N_TD > 0]){
		
		# N_TD_tot <- Run[,N_TD]
		# N_TD_sum <- 0

        # merge Catalogs with row
		Row <- Catalogs[Run][order(as.numeric(gsub(".*[.]([0-9]*)$","\\1",PointSensor)))]
        # TODO: get vDep from Run?
		vdep <- vd[vd_index]

		n <- nrow(Row)
		indCats <- Row[,Cat.Name]
		uniqueCats <- unique(indCats)
		nCats <- length(uniqueCats)
		cat("Sensor", Row[1, Sensor], "/ Source", Row[1, Source],"\n")
		Steps <- unique(round(seq(1,n,length.out=10))) %w/o% n
		N0 <- Row[1,N0]
		cName <- ""
		cSeed <- -1
		Src <- Sources[Sources[,1] %chin% Row[1,Source],]

		Ci <- vector(mode="list",length=n)
		UVW <- vector(mode="list",length=nCats)
		names(UVW) <- uniqueCats

        # data.frame-like or other - i.e. does spatial name column exist?
        if (!is.null(names(vdSpatial[[1]])) && all(names(vdSpatial[[1]]) %in% vdSpatial[[2]][, 1])) {
            # get vDepSpatial column names
            nms_Spatial <- names(vdSpatial[[1]])
            # get corresponding vDep and Sources obj.
            vd_Spatial <- Run[, I(mget(vdSpatial[[1]]))]
        } else {
            # get vDepSpatial column names
            nms_Spatial <- unlist(strsplit(Row[, get(vdSpatial[[1]])], split = ','))
            # get corresponding vDep and Sources obj.
            vd_Spatial <- Run[, I(mget(nms_Spatial))]
        }
        # get corresponding Sources
		Src_Spatial <- lapply(nms_Spatial,function(x,y)y[y[,1] %in% x,],y=vdSpatial[[2]])
        # name them
		names(Src_Spatial) <- names(vd_Spatial) <- nms_Spatial
        # exclude if NA values in vDep
        nms_Spatial <- nms_Spatial[!is.na(unlist(vd_Spatial))]

		for(i in seq_len(n)){

			# if(N_TD_sum == N_TD_tot){
			# 	# browser()
			# 	break
			# }

			# progressbar
			if(i %in% Steps)cat(sprintf("\r[%s%s%s] %1.0f%%",paste(rep(">",round(2*(i-1)/n*10)),collapse=""),"|",paste(rep(".",max(0,19-round(2*(i-1)/n*10))),collapse=""),(i-1)/n*100))
			# {xalt <- matrix(0,2,3)
			# xneu <- gc()
			# while(abs(xalt[2,3]-xneu[2,3])>0){xalt<-xneu;xneu <- gc()}
			# }
			# if(cName!=Row[i,Cat.Name] | cSeed!=Row[i,Subset_seed]){
			if(cName!=Row[i,Cat.Name]){
				cName <- Row[i,Cat.Name]
				cSeed <- Row[i,seed]
				Ctlg <- readCatalog(paste0(C.Path,"/",cName))			
				initializeCatalog(Row[i, ], Catalog = Ctlg)
				uvw <- uvw0(Ctlg)
				Cat.N0 <- attr(Ctlg,"N0")
				rotateCatalog(Ctlg,Row[i,WD])
				if(Cat.N0>N0){
					env <- globalenv()
				    oseed <- env$.Random.seed
				    set.seed(cSeed,kind="L'Ecuyer-CMRG")			
					takeSub <- sample.int(Cat.N0,N0)
			        if (is.null(oseed)) {
			            rm(list = ".Random.seed", envir = env)
			        } else {
			            assign(".Random.seed", value = oseed, envir = env)
			        }
					attCat <- attributes(Ctlg)
					indexNew <- 1:N0
					names(indexNew) <- takeSub
					Ctlg <- Ctlg[Traj_ID %in% takeSub,]
					Ctlg[, ":="(
                        Traj_ID = indexNew[as.character(Traj_ID)]
                        )]
					for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref","uvw0")))setattr(Ctlg,ac,attCat[[ac]])
					uvw <- uvw[takeSub,]
				}
				UVW[[cName]] <- uvw
                # prepare wTD2 once per Ctlg
                if (is_spatial) {
                    Ctlg[, ":="(
                        wTD2 = 2 / wTD
                        )]
                } else {
                    Ctlg[, ":="(
                        wTD2 = 2 / wTD,
                        dep_outside = exp(-vDep * 2 / wTD)
                        )]
                }
			}

			tag_inside(Ctlg,Src,CSnsrs[chmatch(Row[i,PointSensor], CSnsrs[, "Point Sensor Name"]),])
			setkey(Ctlg,Traj_ID)
			# browser()

			if (Ctlg[, any(td_inside)]) {
                # get values
                Ci[[i]] <- ci_fun(Ctlg, nms_Spatial, Src_Spatial, CSnsrs, vdep, Row)
			}
		}
		cat(paste0("\r[",paste0(rep(">",20),collapse=""),"] 100%\n"))

		# rm(Catalogs,CSnsrs,Ctlg,Src,uvw,Run,Row)
		
        .record_now()

		# browser()

		if(n>1){

			# check nulls
			is_null <- sapply(Ci,is.null)
			# # erweitere UVW:
			# is_nUVW <- sapply(UVW,is.null)
			# UVW[is_nUVW] <- UVW[!is_nUVW][1]

			# weights:
			wts <- rep(2,n)
			wts[c(1,n)] <- 1
			rwts <- wts/sum(wts)

			Cmean <- rep(0,n)
			Cmean[!is_null] <- sapply(Ci[!is_null],function(x)x[,sum(CE)])/N0
			# individual CEs
			CEs <- mapply(function(i,x,y){
				out <- rep(0,N0)
				if(i){
					out
				} else {
					out[x[,Traj_ID]] <- x[,CE]
					out - y
				}
			},i=is_null,x=Ci,y=Cmean,SIMPLIFY=TRUE)

			# Us per Cat
			Us <- do.call(cbind,lapply(UVW,function(x){
				x[,"u0"] - mean(x[,"u0"])
			}))
			# Um per Cat
			Um <- do.call(cbind,lapply(UVW,function(x){
				x[,"u0"]
			}))
			# Vs per Cat
			Vs <- do.call(cbind,lapply(UVW,function(x){
				x[,"v0"]
			}))
			# Ws per Cat
			Ws <- do.call(cbind,lapply(UVW,function(x){
				x[,"w0"]
			}))
			rwm <- outer(rwts,rwts)/sqrt(N0)
			Out <- data.frame(
				CE=sum(Cmean*rwts),
				CE_se=sqrt(sum(cov(CEs)*rwm)),
				uCE=sum(colMeans(Us[,indCats[!is_null],drop=FALSE]*CEs[,!is_null,drop=FALSE])*rwts[!is_null]),			
				uCE_se=sqrt(sum(cov(Us[,indCats]*CEs)*rwm)),			
				vCE=sum(colMeans(Vs[,indCats[!is_null],drop=FALSE]*CEs[,!is_null,drop=FALSE])*rwts[!is_null]),
				vCE_se=sqrt(sum(cov(Vs[,indCats]*CEs)*rwm)),
				wCE=sum(colMeans(Ws[,indCats[!is_null],drop=FALSE]*CEs[,!is_null,drop=FALSE])*rwts[!is_null]),
				wCE_se=sqrt(sum(cov(Ws[,indCats]*CEs)*rwm)),
				UCE=sum(colMeans(Um[,indCats[!is_null],drop=FALSE]*t(t(CEs[,!is_null,drop=FALSE]) + Cmean[!is_null]))*rwts[!is_null])
			)

		} else {
			# browser()
			Ci_v <- rep(0,N0)
			Ci_v[Ci[[1]][,Traj_ID]] <- Ci[[1]][,CE]
			Out <- data.frame(
				CE=(CE <- mean(Ci_v)),
				CE_se=sd(Ci_v)/sqrt(N0),
				uCE=mean((UVW[[1]][,"u0"]-mean(UVW[[1]][,"u0"]))*(Ci_v - CE)),
				uCE_se=sd((UVW[[1]][,"u0"]-mean(UVW[[1]][,"u0"]))*(Ci_v - CE))/sqrt(N0),
				vCE=mean(UVW[[1]][,"v0"]*(Ci_v - CE)),
				vCE_se=sd(UVW[[1]][,"v0"]*(Ci_v - CE))/sqrt(N0),
				wCE=mean(UVW[[1]][,"w0"]*(Ci_v - CE)),
				wCE_se=sd(UVW[[1]][,"w0"]*(Ci_v - CE))/sqrt(N0),
				UCE=mean(UVW[[1]][,"u0"]*Ci_v)
			)		
		}
	} else {
		Out <- data.frame(
			CE=0,
			CE_se=NA_real_,
			uCE=0,
			uCE_se=NA_real_,
			vCE=0,
			vCE_se=NA_real_,
			wCE=0,
			wCE_se=NA_real_,
			UCE=0
		)		
	}

    # record gc/mem
    setattr(Out, 'cpu_mem', .record_now(reset = TRUE))

	return(Out)	
}

fill_Ci_homogeneous <- function(Cat, ...) {
    Cat[, dep := 1][(!td_inside), dep := dep_outside]
    # return Traj_ID & CE
    Cat[Traj_ID %in% Traj_ID[(td_inside)],
        .(CE = sum(as.numeric(td_inside) * cumprod(dep) * wTD2))
        , by = Traj_ID]
}

fill_Ci_spatial <- function(Cat, nms_spatial, src_spatial, csnsrs, vd, run_row) {
    # assign for looping over spatial
    Cat[, tagInsideQ := td_inside]
    ### spatially inhomogeneous vdep:
    Cat[, vDep := vd]
    for (j in nms_spatial) {
        tag_inside(Cat, src_spatial[[j]], 
            csnsrs[chmatch(run_row[i, PointSensor], csnsrs[, "Point Sensor Name"]), ]
        )
        Cat[(td_inside), vDep := vd_Spatial[[j]]]
    }
    # TODO: allow for deposition inside source
    #   -> set vdep to 0 inside and remove subsetting below
    #   -> set vdep to 0 before spatial vdep such that it can be set to non
    #           zero through spatial...
    Cat[, dep := 1][(!tagInsideQ), dep := exp(-vDep * wTD2)]
    # return Traj_ID & CE
    Cat[Traj_ID %in% Traj_ID[(tagInsideQ)],
        .(CE = sum(as.numeric(tagInsideQ) * cumprod(dep) * wTD2))
        , by = Traj_ID]
}
