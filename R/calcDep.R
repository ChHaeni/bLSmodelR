.calcDep_Wrapper <- function(RunElement, variables = 'CE', spatial = FALSE) {
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
        variables = variables
    )
    # update vd index
	out[, 'vd_index'] <-  RunElement[, vd_index]
    out
}

.calcDep <- function(Run, Catalogs, C.Path, Sources, CSnsrs, vd, vdSpatial,
    is_spatial, variables = 'CE') {

    # which additional variables?
    which_vars <- c('uCE', 'vCE', 'wCE') %in% variables

    ci_fun <- if (is_spatial) {
        bLSmodelR:::fill_Ci_spatial 
    } else {
        bLSmodelR:::fill_Ci_homogeneous
    }

    # record gc/mem
    .record_now(start = TRUE)

	# on.exit(browser())
	vd_index <- Run[, vd_index]
	# browser()	
	if (Run[, N_TD > 0]) {
		
		# N_TD_tot <- Run[,N_TD]
		# N_TD_sum <- 0

        # merge Catalogs with row
		Row <- Catalogs[Run][order(as.numeric(gsub(".*[.]([0-9]*)$","\\1",PointSensor)))]
        # TODO: get vDep from Run?
		vdep <- vd[vd_index]

		nr <- nrow(Row)
		indCats <- Row[,Cat.Name]
		uniqueCats <- unique(indCats)
		nCats <- length(uniqueCats)
		cat("Sensor", Row[1, Sensor], "/ Source", Row[1, Source],"\n")
		Steps <- unique(round(seq(1,nr,length.out=10))) %w/o% nr
		N0 <- Row[1,N0]
		cName <- ""
		cSeed <- -1
		Src <- Sources[Sources[,1] %chin% Row[1,Source],]

		Ci <- vector(mode="list",length=nr)
		UVW <- vector(mode="list",length=nCats)
		names(UVW) <- uniqueCats

        if (is_spatial) {
            # data.frame-like or other - i.e. does spatial name column exist?
            if (!is.null(names(vdSpatial[[1]])) && 
                    all(names(vdSpatial[[1]]) %in% vdSpatial[[2]][, 1])
                ) {
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
        }

		for(i in seq_len(nr)){

			# if(N_TD_sum == N_TD_tot){
			# 	# browser()
			# 	break
			# }

			# progressbar
			if(i %in% Steps)cat(sprintf("\r[%s%s%s] %1.0f%%",paste(rep(">",round(2*(i-1)/nr*10)),collapse=""),"|",paste(rep(".",max(0,19-round(2*(i-1)/nr*10))),collapse=""),(i-1)/nr*100))
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
                Ci[[i]] <- ci_fun(Ctlg, nms_Spatial, Src_Spatial, CSnsrs, vdep, 
                    vd_Spatial, Row[i, PointSensor])
			}
		}
		cat(paste0("\r[",paste0(rep(">",20),collapse=""),"] 100%\n"))

		# rm(Catalogs,CSnsrs,Ctlg,Src,uvw,Run,Row)
		
        .record_now()

		# browser()

		if (nr > 1) {

			# check nulls
			not_null <- !sapply(Ci, is.null)

			# weights:
			wts <- rep(2, nr)
			wts[c(1, nr)] <- 1
			rwts <- wts / sum(wts)
			orwts <- outer(rwts, rwts) / sqrt(N0)

            # average CE + UCE
            CE_mean <- 0
            UCE_mean <- 0
            for (i_n in which(not_null)) {
                CE_mean <- CE_mean + Ci[[i_n]][, sum(CE)] * rwts[i_n]
                UCE_mean <- UCE_mean + Ci[[i_n]][, 
                    sum(UVW[[indCats[i_n]]][Traj_ID, 'u0'] * CE)] * rwts[i_n]
            }
            CE_mean <- CE_mean / N0
            UCE_mean <- UCE_mean / N0

            # initialize C Matrix
            c_matrix <- matrix(-CE_mean, nrow = N0, ncol = nr)

            # fill C Matrix:
            for (j in which(not_null)) {
                c_matrix[Ci[[j]][, Traj_ID], j] <- Ci[[j]][, CE - CE_mean]
            }

            if (which_vars[1]) {
                # uCE + SE
                Us <- matrix(
                    UVW[[uniqueCats[1]]][, 'u0'] - mean(UVW[[uniqueCats[1]]][, 'u0'])
                    , nrow = N0, ncol = nr)
                for (ic in uniqueCats[-1]) {
                    Us[, indCats == ic] <- UVW[[ic]][, 'u0'] - mean(UVW[[ic]][, 'u0'])
                }
                uvwCE <- c_matrix * Us
                uCE_add <- sum(colSums(uvwCE) * rwts) / N0
                uCE_se_add <- sqrt(sum(cov(uvwCE) * orwts) / N0)
            } else {
                uCE_add <- uCE_se_add <- NA_real_
                if (any(which_vars[2:3])) {
                    Us <- matrix(0.0, nrow = N0, ncol = nr)
                }
            }

            if (which_vars[2]) {
                # vCE + SE
                for (ic in uniqueCats) {
                    Us[, indCats == ic] <- UVW[[ic]][, 'v0']
                }
                uvwCE <- c_matrix * Us
                vCE_add <- sum(colSums(uvwCE) * rwts) / N0
                vCE_se_add <- sqrt(sum(cov(uvwCE) * orwts) / N0)
            } else {
                vCE_add <- vCE_se_add <- NA_real_
            }

            if (which_vars[3]) {
                # uCE + SE
                for (ic in uniqueCats) {
                    Us[, indCats == ic] <- UVW[[ic]][, 'w0']
                }
                uvwCE <- c_matrix * Us
                wCE_add <- sum(colSums(uvwCE) * rwts) / N0
                wCE_se_add <- sqrt(sum(cov(uvwCE) * orwts) / N0)
            } else {
                wCE_add <- wCE_se_add <- NA_real_
            }

			Out <- data.frame(
                # CE
				CE = CE_mean,
				CE_se = sqrt(sum(cov(c_matrix) * orwts) / N0),
                # uCE
				uCE = uCE_add,
				uCE_se = uCE_se_add,
                # vCE
				vCE = vCE_add,
				vCE_se = vCE_se_add,
                # wCE
				wCE = wCE_add,
				wCE_se = wCE_se_add,
                # UCE
				UCE = UCE_mean
			)

		} else {
			Ci_v <- rep(0, N0)
			Ci_v[Ci[[1]][, Traj_ID]] <- Ci[[1]][, CE]
            CE_mean <- mean(Ci_v)
            dCi <- Ci_v - CE
            dui <- UVW[[1]][, 'u0'] - mean(UVW[[1]][, 'u0'])
			Out <- data.frame(
				CE = CE_mean,
				CE_se = sd(Ci_v) / sqrt(N0), 
				uCE = mean(dui * dCi), 
				uCE_se = sd(dui * dCi) / sqrt(N0), 
				vCE = mean(UVW[[1]][, "v0"] * dCi), 
				vCE_se = sd(UVW[[1]][, "v0"] * dCi) / sqrt(N0), 
				wCE = mean(UVW[[1]][, "w0"] * dCi), 
				wCE_se = sd(UVW[[1]][, "w0"] * dCi) / sqrt(N0), 
				UCE = mean(UVW[[1]][, "u0"] * Ci_v)
			)		
		}
	} else {
		Out <- data.frame(
			CE = 0, 
			CE_se = NA_real_, 
			uCE = 0, 
			uCE_se = NA_real_, 
			vCE = 0, 
			vCE_se = NA_real_, 
			wCE = 0, 
			wCE_se = NA_real_, 
			UCE = 0
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

fill_Ci_spatial <- function(Cat, nms_spatial, src_spatial, csnsrs, vd, vd_spatial, 
    point_sensor) {
    # assign for looping over spatial
    Cat[, tagInsideQ := td_inside]
    ### spatially inhomogeneous vdep:
    Cat[, vDep := vd]
    for (j in nms_spatial) {
        tag_inside(Cat, src_spatial[[j]], 
            csnsrs[chmatch(point_sensor, csnsrs[, "Point Sensor Name"]), ]
        )
        Cat[(td_inside), vDep := vd_spatial[[j]]]
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
