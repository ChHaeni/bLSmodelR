rebuildCatListFile <- function(C.Path, fromScratch = FALSE) {
    # catalog flag -> same as readCatalog/writeCatalog
    catalog_flag <- 'A'
    # check if cats exist & list them
	Existing <- dir(C.Path, pattern = "Cat_Zm.*_[0-9]{14}$")
    # define path to CatList file
	Catfile <- paste0(C.Path, "/.cats")
    # check length of Existing
	if (length(Existing) > 0) {
        # full path names
        ExistingFull <- paste(C.Path, Existing, sep = "/")
        # get CatList (read file)
		if (file.exists(Catfile)) {
            # try to read file
            CatList <- try(qread(Catfile, strict = TRUE), silent = TRUE)
            # try again on error
            tic <- Sys.time()
            while (inherits(CatList, 'try-error') && as.numeric(Sys.time() - tic, units = 'secs') < 20) {
                CatList <- try(qread(Catfile, strict = TRUE), silent = TRUE)
            }
            # improved error message upon failure
            if (inherits(CatList, 'try-error')) {
                stop('rebuildCatListFile: reading file ', Catfile, ' fails!')
            }
            # check catalog flag
            cat_flag <- attr(CatList, 'cat_flag')
            if (is.null(cat_flag) || cat_flag != catalog_flag) {
                CatList <- as.data.table(c(list(a = character(0), b = as.POSIXct(character(0), tz = 'GMT')), rep(list(a = numeric(0)), 13)), stringsAsFactors = FALSE)
                setnames(CatList, c("Name", "mtime", "N0", "ZSens", "Ustar", "L", "Zo", "Su_Ustar", "Sv_Ustar", "bw", "C0", "kv", "A", "alpha", "MaxFetch"))
            } else {
                # check erroneous
                CatList <- CatList[grepl("Cat_Zm.*_[0-9]{14}$", Name)]
                # remove duplicates
                CatList <- unique(CatList, by = 'Name')
            }
		} else {
			CatList <- as.data.table(c(list(a = character(0), b = as.POSIXct(character(0), tz = 'GMT')), rep(list(a = numeric(0)), 13)), stringsAsFactors = FALSE)
			setnames(CatList, c("Name", "mtime", "N0", "ZSens", "Ustar", "L", "Zo", "Su_Ustar", "Sv_Ustar", "bw", "C0", "kv", "A", "alpha", "MaxFetch"))
        }
        # exclude any non-existing
        CatList <- CatList[Name %chin% Existing]
        # set Name as key
        setkey(CatList, Name)
        # exclude modified
        CatList <- CatList[mtime == file.mtime(file.path(C.Path, Name))]
        # check catalogs not in CatList
		if (any(checkCat <- CatList[, !(Existing %chin% Name)])){
            cat('-> updating catalog db...\n')
            # create CatAdd to append at bottom
            nr <- sum(checkCat)
			CatAdd <- setNames(
                as.data.frame(c(list(a = character(nr), b = as.POSIXct(rep(Sys.time(), nr), tz = 'GMT')), rep(list(a = numeric(nr)), 13)), stringsAsFactors = FALSE),
			    c("Name", "mtime", "N0", "ZSens", "Ustar", "L", "Zo", "Su_Ustar", "Sv_Ustar", "bw", "C0", "kv", "A", "alpha", "MaxFetch"))
            # get check index
            check_index <- which(checkCat)
			for(j in seq_along(check_index)){
                cat('\r', j, '/', length(check_index))
                # get i
                i <- check_index[j]
                # read catalog
				CatHeader <- readCatalog(ExistingFull[i], header_only = TRUE)
                # convert from old serialization?
                if (is.null(CatHeader)) {
                    # old qs format
                    Cat <- try(qs::qread(ExistingFull[i], strict = TRUE), silent = TRUE)
                    if (inherits(Cat, 'try-error')) {
                        # old rds format
                        Cat <- try(readRDS(ExistingFull[i]))
                    }
                    if (inherits(Cat, 'try-error')) {
                        # file corrupt
                        file.remove(ExistingFull[i])
                    } else {
                        # save with new binary format
                        writeCatalog(Cat, ExistingFull[i])
                        CatHeader <- attr(Cat, 'header')
                    }
                }
                if (!is.null(CatHeader)) {
                    Head <- unlist(strsplit(CatHeader, "\n"))[-1]
                    Whead <- matrix(as.numeric(gsub(".*[=] ", "", Head)), nrow=1)
                    CatAdd[j, -(1:2)] <- Whead
                    # get file modified
                    CatAdd[j, 2] <- file.mtime(ExistingFull[i])
                    # get file name
                    CatAdd[j, 1] <- Existing[i]
                }
			}
            # bind together
			CatList <- na.omit(rbind(CatList, CatAdd))
            cat('\r                                               \r   done\n')
		}
        # add catalog flag as attribute
        setattr(CatList, 'cat_flag', catalog_flag)
        # try to write - if error retry for 20 seconds
        try_write <- try(qsave(CatList, Catfile, preset = 'uncompressed'), silent = TRUE)
        # loop on error
        time_now <- Sys.time()
        while(inherits(try_write, 'try-error')) {
            # check 20 secs
            if (as.numeric(Sys.time() - time_now, units = 'secs') >= 20) {
                stop("rebuildCatListFile: can't write file ", Catfile,
                    "\nerror message: ", conditionMessage(attr(try_write, 'condition')))
            }
            # wait to continue...
            Sys.sleep(1)
            # try again
            try_write <- try(qsave(CatList, Catfile, preset = 'uncompressed'), silent = TRUE)
        }
	} else {
        ## no catalogs exist
		CatList <- as.data.frame(c(list(a=character(0)), rep(list(a=numeric(0)), 13)), stringsAsFactors=FALSE)
		colnames(CatList) <- c("Name", "N0", "ZSens", "Ustar", "L", "Zo", "Su_Ustar", "Sv_Ustar", "bw", "C0", "kv", "A", "alpha", "MaxFetch")
        # remove file if it exists
		if (file.exists(Catfile)) file.remove(Catfile)
	}
    # return invisible
	invisible(CatList)
}

