deposition <- function(x, vDep, rn = NULL, Sensor = NULL, Source = NULL, 
    vDepSpatial = NULL, ncores = 1, memory_limit = NULL) {

	# argument names better!
	RNG <- RNGkind(kind="L'Ecuyer-CMRG")

	on.exit(
		{
			RNGkind(kind=RNG[1])
		}
	)

  # convert old versions 
  sx <- as.character(substitute(x))
	x <- copy(x)
	setDT(x)
	switchNames(x)
  if(is.null(attr(x, "Version"))){
		warning(paste0("Object '", sx[min(length(sx), 2)], "' has not yet been converted to version 4.2+"))
		convert(x)
  }
  
	ModelInput <- attr(x,"ModelInput")
	Catalogs <- attr(x,"Catalogs")
	Cat.Path <- attr(x,"CatPath")

	if(is.null(rn)){
		Selrn <- x[,unique(rn)]
	} else {
		Selrn <- rn
	}
	if(is.null(Sensor)){
		SelSensor <- x[,unique(Sensor)]
	} else {
		SelSensor <- Sensor
	}
	if(is.null(Source)){
		SelSource <- x[,unique(Source)]
	} else {
		SelSource <- Source
	}
	Run <- x[rn %in% Selrn & Sensor %chin% SelSensor & Source %chin% SelSource,]

	N <- nrow(Run)
	Run[, vd_index := 1:.N]

    # TODO: attach vDep to Run in any case?
	if(is.character(vDep)){
		vDep <- Run[,vDep,with=FALSE][[1]]
	} else {
		vDep <- rep(vDep,N)[seq_len(N)]
	}

	# vDepSpatial
	if(vdSpat <- !is.null(vDepSpatial)){
        # check 'Sources' obj
		if(!inherits(vDepSpatial[[2]],"Sources")){
			stop("Second list element of argument 'vDepSpatial' must be of class 'Sources'!")
		}
        # copy first obj
        vds1 <- copy(vDepSpatial[[1]])
        # get names
		nms <- names(vds1)
        # check if list, data.frame-like or character-vector
        if (is.character(vds1)) {
            # column names (old option, but allow for NA -> new option 2.)
            # check if names are null
            if (is.null(nms)) {
                # assign values as names if missing
                names(vDepSpatial[[1]]) <- nms <- 
                    names(vds1) <- vds1 
            }
            # check if column names exist
            if (!all(vds1 %in% names(Run))) {
                stop('columns ', paste(vds1[!(vds1 %in% names(Run))],collapse=", ")," are missing!")
            }
            # check names:
            if(!all(nms %in% unique(vDepSpatial[[2]][,1]))){
                stop(paste(nms[!(nms %in% unique(vDepSpatial[[2]][,1]))],collapse=", "),": area not defined!")
            }
        } else if (inherits(vds1, 'data.frame')) {
            # new option 1.
            # subset by rn, Sensor & Source if nrow != N
            if (nrow(vds1) != N) {
                # must contain columns to subset
                if (!all(c('rn', 'Sensor', 'Source') %in% nms)) {
                    stop("First list element of argument 'vDepSpatial' has non-matching number of rows.",
                        "\n  Provide columns 'rn', 'Sensor' and 'Source' for subsetting.")
                }
                # subset
                vds1 <- vds1[rn %in% Selrn & Sensor %chin% SelSensor & Source %chin% SelSource, ]
            }
            # remove subsetting column names and change to lower-case
            nms_resid <- setdiff(nms, c('rn', 'Source', 'Sensor'))
            # find numeric columns:
            nm_numeric <- vds1[, which(unlist(lapply(.SD, is.numeric))), .SDcols = nms_resid]
            # find 'Spatial' column: 
            # try if there's only one column containing characters
            nmi_spatial <- vds1[, which(unlist(lapply(.SD, is.character))), .SDcols = nms_resid]
            if (length(nmi_spatial) > 1) {
                # try to be smart
                for (nm in nmi_spatial) {
                    # check if any numeric name can be found
                    if (vds1[, sum(grepl(paste(nms_resid[nm_numeric], collapse = '|'), get(nms_resid[nm]))) == 0]) {
                        # remove from nmi_spatial if not any of the numeric columns can be found
                        nmi_spatial <- setdiff(nmi_spatial, nm)
                    }
                }
                # still not unique?
                if (length(nmi_spatial) > 1) {
                    # tolower -> grep spatial, zones, source, select, regions -> check character -> find sources
                    nmi_spatial <- nmi_spatial[grep('spatial|zones|source|select|regions|areas|patches', tolower(nms_resid[nmi_spatial]))]
                }
            }
            # did we find (only) one spatial column?
            if (length(nmi_spatial) == 0) {
                # none
                stop("Cannot find a column with Source Name entries in the first list element of argument 'vDepSpatial'")
            } else if (length(nmi_spatial) > 1) {
                # not unique
                stop("Cannot find a column with Source Name entries in the first list element of argument 'vDepSpatial'")
            }
            # get spatial column name
            nm_spatial <- nms_resid[nmi_spatial]
            # get sources (split , separated)
            nms_sources <- vds1[, unique(unlist(strsplit(unique(get(nm_spatial)), split = ',')))]
            # don't allow NA values
            if (anyNA(nms_sources)) {
                stop("NA values in Source Name entries (first list element of argument 'vDepSpatial'.",
                    " Use empty string '' if you don't want to include any vDepSpatial patches for individual intervals.")
            }
            # check if columns exist
            if (!all(chk_nms <- nms_sources %in% nms_resid)) {
                stop("Cannot find column(s) ", paste(nms_sources[!chk_nms], collapse = ', '), " in first list element of argument 'vDepSpatial'")
            }
            # rename to paste0('vDepSpatial.', *)
            # spatial column
            for (nm_src in nms_sources) {
                # update spatial source name
                vds1[, (nm_spatial) := sub(paste0('\\<', nm_src, '\\>'), paste0('vDepSpatial.', nm_src), get(nm_spatial))]
            }
            # column names
            setnames(vds1, nms_resid, paste0('vDepSpatial.', nms_resid))
            # attach or merge to Run
            if (all(c('rn', 'Sensor', 'Source') %in% nms)) {
                # merge
                Run <- merge(Run, vds1, by = c('rn', 'Sensor', 'Source'))
            } else {
                # blindly cbind
                Run <- cbind(Run, vds1[, setdiff(names(vds1), c('rn', 'Sensor', 'Source')), with = FALSE])
            }
            # update vDepSpatial[[1]] to contain the column name of 'spatial column' -> get column names from 'spatial column' entries
            vDepSpatial[[1]] <- paste0('vDepSpatial.', nm_spatial)
            # update vDepSpatial[[2]] names
            vDepSpatial[[2]][, 1] <- paste0('vDepSpatial.', vDepSpatial[[2]][, 1])
        } else if (is.list(vds1)) {
            # old option
            # check names:
            if (is.null(nms)) {
                stop("First list element of argument 'vDepSpatial' must be named")
            }
            if(!all(nms %in% unique(vDepSpatial[[2]][,1]))){
                stop(paste(nms[!(nms %in% unique(vDepSpatial[[2]][,1]))],collapse=", "),": area not defined!")
            }
            # extend/recycle
            for(i in nms){
                if(is.character(vds1[[i]])){
                    vds1[[i]] <- Run[, vds1[[i]], with = FALSE][[1]]
                } else {
                    n_rep <- ceiling(N / length(vds1[[i]]))
                    vds1[[i]] <- rep(vds1[[i]], n_rep)[seq_len(N)]
                }
            }
            # attach to Run -> cbind
            Run <- cbind(Run, as.data.table(vds1))
            # replace vDepSpatial[[1]]
            vDepSpatial[[1]] <- setNames(nms, nms)
        } else {
            stop("First list element of argument 'vDepSpatial' must be either a character vector,",
                " a list or an object that inherits from 'data.frame'!")
        }
	}

	### procSensors:
	pSens <- procSensors(ModelInput[["Sensors"]])

	SourceNames <- unique(ModelInput[["Sources"]][,1])
	# fracDepList <- vector(length(SourceNames),mode="list")
	# names(fracDepList) <- SourceNames
	# vDepList <- fracDepList
	# if(length(fracDepInside)==1){
	# 	fracDepList[] <- fracDepInside
	# } else {
	# 	if(!all(SourceNames%in%names(fracDepInside)))stop("Missing sources names in argument fracDepInside!\nSource Names are: \"",paste(SourceNames,collapse="\", \""),"\"")
	# 	fracDepList[] <- fracDepInside[SourceNames]
	# }
	# if(length(vDepInside)==1){
	# 	vDepList[] <- vDepInside
	# } else {
	# 	if(!all(SourceNames%in%names(vDepInside)))stop("Missing sources names in argument vDepInside!\nSource Names are: \"",paste(SourceNames,collapse="\", \""),"\"")
	# 	vDepList[] <- vDepInside[SourceNames]
	# }
	# vDepList <- as.data.frame(vDepList)

	setkey(Run,rn,Sensor)

	if (inherits(ncores, 'cluster') || (ncores != 1 && N > 1)) {


        if (inherits(ncores, 'cluster')) {
            cl <- ncores
            ncores <- length(cl)
        } else {
			if (ncores < 1) ncores <- parallel::detectCores(TRUE, FALSE)
			on.exit(parallel::stopCluster(cl), add = TRUE)
            cl <- .makePSOCKcluster(ncores, memory_limit = memory_limit)
			data.table::setDTthreads(ncores)
		}
        parallel::clusterEvalQ(cl, data.table::setDTthreads(1L))

        # get N_TD > 0
        ntd_g0 <- Run[, which(N_TD > 0)]

        if (length(ntd_g0) > 0) {

            # sort by N_TD
            ntd_order <- Run[ntd_g0, order(N_TD)]

            ### try to qsave/qread, instead of direct traffic?
            # get temporary file name
            # NOTE: what if tempdir is not accessable (e.g. PRIME?)
            t_file <- tempfile(pattern = 'parbls_input_')

            # distribute to clusters
            InputFiles <- lapply(seq_along(cl), function(x, run) {
                    # file name
                    f_name <- paste0(t_file, '_', x)
                    # qsave to tempfile
                    qs::qsave(list(
                            RunList = run[seq.int(from = x, to = nrow(run), by = length(cl)), ],
                            Catalogs = Catalogs, 
                            Cat.Path = Cat.Path, 
                            Sources = ModelInput[['Sources']], 
                            Sensors = pSens$'Calc.Sensors', 
                            vDep = vDep, 
                            vDepSpatial = vDepSpatial), 
                        file = f_name,
                        preset = 'uncompressed'
                    )
                    # return file name
                    f_name
                }, run = as.data.frame(Run[ntd_g0[ntd_order], ]))

            # run parallel
            cat("Parallel computing deposition corrected C/E...\nThis will take a moment...\n")
            if(vdSpat){
                ResFiles <- try(parallel::clusterApply(cl, InputFiles, .calcDep_Wrapper, spatial = TRUE), silent = TRUE)
            } else {
                ResFiles <- try(parallel::clusterApply(cl, InputFiles, .calcDep_Wrapper), silent = TRUE)
            }
            # check try-error
            if (inherits(ResFiles, 'try-error')) {
                parallel::stopCluster(cl)
                stop('parallel computing returned the following error message:\n',
                    attr(ResFiles, 'condition')[['message']])
            }

            # read files
            OutList <- lapply(ResFiles, qs::qread)
            # clean up
            file.remove(unlist(ResFiles))

            # bind together (we do not need to order because of merge)
            Out <- rbind(
                # N_TD > 0 from clusters
                rbindlist(OutList),
                # N_TD == 0
                Run[N_TD == 0, .(CE, CE_se, uCE, uCE_se, vCE, vCE_se, wCE, wCE_se, UCE, vd_index)]
            )

        } else {

            # all N_TD == 0
            Out <- Run[, .(CE, CE_se, uCE, uCE_se, vCE, vCE_se, wCE, wCE_se, UCE, vd_index)]

        }

	} else {

		OutList <- vector(N,mode="list")
		if(vdSpat){
			for(i in seq_len(N)){
				cat(i,"/",N,":\n")
				OutList[[i]] <- .calcDep_Spatial(Run[i,],Catalogs,Cat.Path,ModelInput[["Sources"]],pSens$"Calc.Sensors",vDep,vDepSpatial)
			}
		} else {
			for(i in seq_len(N)){
				cat(i,"/",N,":\n")
				OutList[[i]] <- .calcDep(Run[i,],Catalogs,Cat.Path,ModelInput[["Sources"]],pSens$"Calc.Sensors",vDep,vDepSpatial)
			}			
		}
		Out <- rbindlist(OutList)[, vd_index := seq_len(N)]
	}

	setnames(Out,paste0(names(Out),"_Dep"))
	Out <- merge(Run, Out, by.x = "vd_index", by.y = "vd_index_Dep")
	Out[, vd_index := NULL]

    # pass old attributes
    # TODO: add function bls_attributes -> CalcSteps, CatPath, Catalogs,
    #           ModelInput, Version, ModelRunTime, vDep
    for (att in c('CalcSteps', 'CatPath', 'Catalogs', 'ModelInput', 'Version', 'ModelRunTime')) {
        setattr(Out, att, attr(Run, att))
    }

    # add new attributes
	setattr(Out,"vDep",list(vDep=vDep,vDepSpatial=vDepSpatial))
	setattr(Out,"class",c("deposition",class(Out)))

	return(Out)
}
