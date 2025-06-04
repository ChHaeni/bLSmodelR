.onLoad <- function(libname, pkgname) {
	 invisible(
         reg.finalizer(
        e = parent.env(environment()),
        f = function(env){
            eval(env, {
                bLSmodelR::cleanTemporary()
            })
        },
        onexit = TRUE))
}
.onAttach <- function(libname, pkgname){
    packageStartupMessage("\n#################################")
    packageStartupMessage(paste0(" This is bLSmodelR version ",packageVersion("bLSmodelR")))
    packageStartupMessage(paste0(" last updated ",packageDescription("bLSmodelR")[["Date"]]))
    packageStartupMessage("#################################\n")
    # set some options
    options(
        # parent directory of job dir
        bls.slurm.jobdir = getOption('bls.slurm.jobdir', file.path(Sys.getenv('HOME'), '.slurm')),
        # exclude partitions?
        bls.slurm.exclude.partition = getOption('bls.slurm.exclude.partition', '')
    )
}
.inifu <- function(n,theta){
	# http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
	Id <- diag(n)
	x1    <- rnorm(n)        # fixed given data
	x1 <- (x1 - mean(x1))#/sd(x1)
	x2    <- rnorm(n)      # new random data
	Xctr   <- cbind(x1, x2 - mean(x2))
	Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
	P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
	x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
	Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
	Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
	x3    <- rnorm(n)      # new random data
	Xctr     <- cbind(Y, x3 - mean(x3))         # matrix
	Q    <- qr.Q(qr(Xctr[ , 1:2, drop=FALSE]))      # QR-decomposition, just matrix Q
	P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
	x3o  <- (Id-P) %*% Xctr[ , 3]                 # x2ctr made orthogonal to x1ctr
	Xc2  <- cbind(Xctr[ , 1], x3o)                # bind to matrix
	Y2    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
	x2 <- Y[ , 2]    # final new vector
	x3 <- Y2[ , 2] + (1 / tan(theta)) * Y[ , 1]    # final new vector
	return(cbind(x1/sd(x1),x2/sd(x2),x3/sd(x3)))
}

.CheckCatMatches <- function(int_ext, cat_list, tol, tol_lower, tol_upper){
    # rename
    setnames(cat_list, names(cat_list)[-1], paste0("Cat_", names(cat_list)[-1]))
    # create Cat_Sensor_Swustar
    cat_list[, Cat_Sensor_Swustar := round(calcsigmaW(1, Cat_ZSens / Cat_L, Cat_bw), 3)]

    # check matching catalogs 
    Key <- int_ext[, {
        cat('\r\r** check grouped intervals:', .GRP, '/', .NGRP)
        # check match within Tolerances & MaxFetch
        cat_list[
            Cat_MaxFetch >= .BY[['MaxFetch']]
        ][
            Cat_ZSens >= .BY[['z_lo']] & Cat_ZSens <= .BY[['z_up']]
        ][
            Cat_Sensor_Swustar >= .BY[['sw_lo']] & Cat_Sensor_Swustar <= .BY[['sw_up']]
        ][
            Cat_L >= .BY[['l_lo']] & Cat_L <= .BY[['l_up']]
        ][
            Cat_Zo >= .BY[['z0_lo']] & Cat_Zo <= .BY[['z0_up']]
        ][
            Cat_Su_Ustar >= .BY[['sUu_lo']] & Cat_Su_Ustar <= .BY[['sUu_up']]
        ][
            Cat_Sv_Ustar >= .BY[['sVu_lo']] & Cat_Sv_Ustar <= .BY[['sVu_up']]
        ][, .(
            index = paste(row, collapse = '-'),
            rows = list(row),
            cat_name = Name,
            cat_n0 = Cat_N0,
            Cat_ZSens = Cat_ZSens, 
            Cat_L = Cat_L,
            Cat_Zo = Cat_Zo,
            Cat_Su_Ustar = Cat_Su_Ustar,
            Cat_Sv_Ustar = Cat_Sv_Ustar,
            Cat_Sensor_Swustar = Cat_Sensor_Swustar,
            Cat_bw = Cat_bw,
            Cat_MaxFetch = Cat_MaxFetch,
            Cat_C0 = Cat_C0,
            Cat_alpha = Cat_alpha,
            Cat_A = Cat_A,
            Cat_kv = Cat_kv,
            devZSens = abs(Cat_ZSens / SensorHeight[1] - 1) / tol[1], 
            devL = abs(Cat_L / L[1] - 1) / tol[2], 
            devZo = abs(Cat_Zo / Zo[1] - 1) / tol[3], 
            devsUu = abs(Cat_Su_Ustar / sUu[1] - 1) / tol[4], 
            devsVu = abs(Cat_Sv_Ustar / sVu[1] - 1) / tol[5], 
            devSensor_Swustar = abs(Cat_Sensor_Swustar / Sensor_Swustar[1] - 1) / tol[6], 
            devN0 = Cat_N0 - .BY[['N0']]
        )]
    }, by = 
        .(
            kv, A, alpha, MaxFetch, N0,
            z_lo = SensorHeight * tol_lower[, 'Sensor Height'],
            z_up = SensorHeight * tol_upper[, 'Sensor Height'],
            l_lo = ifelse(L < 0, L * tol_upper[, 'L'], L * tol_lower[, 'L']),
            l_up = ifelse(L < 0, L * tol_lower[, 'L'], L * tol_upper[, 'L']),
            z0_lo = Zo * tol_lower[, 'Zo'],
            z0_up = Zo * tol_upper[, 'Zo'],
            sw_lo = Sensor_Swustar * tol_lower[, 'SigmaW/Ustar'],
            sw_up = Sensor_Swustar * tol_upper[, 'SigmaW/Ustar'],
            sUu_lo = sUu * tol_lower[, 'SigmaU/Ustar'],
            sUu_up = sUu * tol_upper[, 'SigmaU/Ustar'],
            sVu_lo = sVu * tol_lower[, 'SigmaV/Ustar'],
            sVu_up = sVu * tol_upper[, 'SigmaV/Ustar']
        )
    ][
        devZSens <= 1.0000001 &
        devL <= 1.0000001 &
        devZo <= 1.0000001 &
        devsUu <= 1.0000001 &
        devsVu <= 1.0000001 &
        devSensor_Swustar <= 1.0000001, sumDev := devZSens + devL + devZo + devsUu + devsVu + devSensor_Swustar
    ]
    cat('\n')

    if(nrow(Key) > 0){

        # add for cat
        int_ext[, Cat.extend := FALSE]

        # note to myself: apply biases in sumDev someday (see below for explanation)
        Key[devN0 > 0, sumDev := sumDev + 0.001]
        Key[devN0 < 0, sumDev := sumDev + 6 + log(-devN0)]

        # get best match
        Key[, {
            cat('\r\r** check best matches:', .GRP, '/', .NGRP)
            # which min dev
            ind <- which.min(sumDev)
            # assign catalog name & set cat.exists & cat.calc
            int_ext[row %in% rows[[ind]], ':='(
                    Cat.exists = TRUE,
                    Calc.N0 = max(N0),
                    Cat.calc = c(cat_n0[ind] < max(N0), rep(FALSE, .N - 1)), 
                    Cat.Name = cat_name[ind],
                    Cat.extend = c(cat_n0[ind] < max(N0), rep(FALSE, .N - 1)),
                    Calc.ZSens = Cat_ZSens[ind], 
                    Calc.L = Cat_L[ind],
                    Calc.Zo = Cat_Zo[ind],
                    Calc.Su_Ustar = Cat_Su_Ustar[ind],
                    Calc.Sv_Ustar = Cat_Sv_Ustar[ind],
                    Calc.Sensor_Swustar = Cat_Sensor_Swustar[ind],
                    Calc.bw = Cat_bw[ind],
                    Calc.MaxFetch = Cat_MaxFetch[ind],
                    Calc.C0 = Cat_C0[ind],
                    Calc.alpha = Cat_alpha[ind],
                    Calc.A = Cat_A[ind],
                    Calc.kv = Cat_kv[ind]
                )]
            NULL
        }, by = index]
        cat('\n')

    }

    invisible(int_ext)

}

.getKey <- function(int_ext, row_list, tol, tol_lower, tol_upper){
    # subset only missing
    int_miss <- int_ext[!(Cat.exists)]
    # add limits
    int_miss[, ':='(
        z_lo = SensorHeight * tol_lower[, 'Sensor Height'],
        z_up = SensorHeight * tol_upper[, 'Sensor Height'],
        l_lo = ifelse(L < 0, L * tol_upper[, 'L'], L * tol_lower[, 'L']),
        l_up = ifelse(L < 0, L * tol_lower[, 'L'], L * tol_upper[, 'L']),
        z0_lo = Zo * tol_lower[, 'Zo'],
        z0_up = Zo * tol_upper[, 'Zo'],
        sw_lo = Sensor_Swustar * tol_lower[, 'SigmaW/Ustar'],
        sw_up = Sensor_Swustar * tol_upper[, 'SigmaW/Ustar'],
        sUu_lo = sUu * tol_lower[, 'SigmaU/Ustar'],
        sUu_up = sUu * tol_upper[, 'SigmaU/Ustar'],
        sVu_lo = sVu * tol_lower[, 'SigmaV/Ustar'],
        sVu_up = sVu * tol_upper[, 'SigmaV/Ustar']
    )]
    # prepare Key
    Key <- int_miss[, {
        cat('\r\r** check grouped intervals:', .GRP, '/', .NGRP)
        # check matches within Tolerances & MaxFetch
        # -> .SD provides TD catalog, row_list provides possible matching rows
        row_list[
            Cat_MaxFetch <= .BY[['MaxFetch']]
        ][
            Cat_ZSens >= .BY[['z_lo']] & Cat_ZSens <= .BY[['z_up']]
        ][
            Cat_Sensor_Swustar >= .BY[['sw_lo']] & Cat_Sensor_Swustar <= .BY[['sw_up']]
        ][
            Cat_L >= .BY[['l_lo']] & Cat_L <= .BY[['l_up']]
        ][
            Cat_Zo >= .BY[['z0_lo']] & Cat_Zo <= .BY[['z0_up']]
        ][
            Cat_Su_Ustar >= .BY[['sUu_lo']] & Cat_Su_Ustar <= .BY[['sUu_up']]
        ][
            Cat_Sv_Ustar >= .BY[['sVu_lo']] & Cat_Sv_Ustar <= .BY[['sVu_up']]
        ][, {
            # row: -> identical intervals which can read the same catalog
            # Zeile: -> index for possible catalog rows for given row(s)
            # first_row: provides TD catalog
            .(
                # index = paste(row, collapse = '-'),
                first_row = row[1],
                matched_rows = list(Zeile),
                n_matched = .N,
                cat_name = Cat.Name[1],
                devZSens = max(abs(Cat_ZSens / SensorHeight[1] - 1) / tol[1]), 
                devL = max(abs(Cat_L / L[1] - 1) / tol[2]), 
                devZo = max(abs(Cat_Zo / Zo[1] - 1) / tol[3]), 
                devsUu = max(abs(Cat_Su_Ustar / sUu[1] - 1) / tol[4]), 
                devsVu = max(abs(Cat_Sv_Ustar / sVu[1] - 1) / tol[5]), 
                devSensor_Swustar = max(abs(Cat_Sensor_Swustar / Sensor_Swustar[1] - 1) / tol[6]), 
                devN0 = max(Cat_N0 - .BY[['N0']])
            )
        }]
    }, by = .(kv, A, alpha, MaxFetch, N0, 
        z_lo, z_up, l_lo, l_up, z0_lo, z0_up, 
        sw_lo, sw_up, sUu_lo, sUu_up, sVu_lo, sVu_up)][
                    devZSens <= 1.0000001 &
                    devL <= 1.0000001 &
                    devZo <= 1.0000001 &
                    devsUu <= 1.0000001 &
                    devsVu <= 1.0000001 &
                    devSensor_Swustar <= 1.0000001, 
            sumDev := devZSens + devL + devZo + devsUu + devsVu + devSensor_Swustar
    ]
    cat('\n')
    Key
}

.CheckCrossMatches <- function(int_ext, Key){
    # check cross-matches
    if (nrow(Key) == 0) {
        int_ext[!(Cat.exists), Cat.calc := TRUE]
    } else {
        # apply bias to catalogs with fewer trajectories, 
        #   and make exact N0 matches (just) preferable
        Key[devN0 > 0, sumDev := sumDev + 0.001]
        Key[devN0 < 0, sumDev := sumDev + 6 + log(-devN0)]
        # sort descending
        setorder(Key, -n_matched, sumDev)
        # add helper column (don't count already checked...)
        int_ext[, cat_calc := Cat.exists]
        # loop is sorted (see ?data.table -> by)
        Key[, {
            cat('\r\r** check best cross-matches:', .GRP, '/', .NGRP)
            # check if not yet calculated
            if (int_ext[.BY[[1]], !cat_calc]) {
                # get cross-matching rows which are not yet covered
                icm <- int_ext[matched_rows[[1]]][!(cat_calc), row]
                # set cat_calc to TRUE (to skip in next checks)
                int_ext[icm, cat_calc := TRUE]
                # assign all rows except first_row to FALSE and fix catalog names
                int_ext[icm, 
                    c('Cat.calc', 'Cat.Name') := .(FALSE, cat_name)]
                # fix first row (Cat.calc must be TRUE)
                int_ext[.BY[[1]], Cat.calc := TRUE]
            }
            NULL
        }, by = first_row]
        cat('\n')
    }
    invisible(int_ext)
}

.MaxFetchWrapper <- function(Int_Ext, p_Sens, Input_List){
	Int_Ext[MaxFetch < 0,
		MaxFetch := {
			nmSens <- unlist(strsplit(Sensor,split=","))
			nmSous <- unlist(strsplit(Source,split=","))
			Sens <- structure(p_Sens$Calc.Sensors[p_Sens$Calc.Sensors[, "Point Sensor Name"] %in% nmSens,
				c("x-Coord (m)", "y-Coord (m)")],class="data.frame")
			Sous <- structure(Input_List[["Sources"]][Input_List[["Sources"]][,1] %in% nmSous,2:3],class="data.frame")
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

### add helpers to log memory usage on slurm/parallel workers
# .get_object_sizes <- function(envir = parent.frame()) {
#     objs <- ls(envir = envir)
#     obj.sizes <- lapply(objs, function(x) {
#         object.size(get(x, envir = envir))
#                 })
#     objs_size <- unlist(lapply(obj.sizes, format,
#             units = 'auto', standard = 'SI'))
#     total_size <- format(structure(
#             sum(unlist(obj.sizes)), class = 'object_size'),
#             units = 'auto', standard = 'SI')
#     structure(data.frame(
#         obj = c(objs, '===', 'total:'),
#         size = c(objs_size, '===', total_size)
#     ), total = sum(unlist(obj.sizes)))
# }
# .record_object_sizes <- function(start = FALSE, reset = FALSE, envir = parent.frame()) {
#     if (start) {
#         otab_old <- structure(0, total = 0)
#     } else {
#         otab_old <- getOption('.bls_obj_sizes', 
#             default = structure(0, total = 0))
#     }
#     obj_table <- .get_object_sizes(envir)
#     if (attr(otab_old, 'total') > attr(obj_table, 'total')) {
#         obj_table <- otab_old
#     }
#     if (reset) {
#         options('.bls_obj_sizes' = NULL)
#     } else {
#         options('.bls_obj_sizes' = obj_table)
#     }
#     invisible(obj_table)
# }
# .record_gc_mem <- function(start = FALSE, reset = FALSE) {
#     if (start) options('.bls_gc_mem' = NULL)
#     new_gc <- gc()
#     old_gc <- getOption('.bls_gc_mem', matrix(0, nrow = 2, ncol = 6))
#     for (i in c(2, 4, 6)) {
#         if (new_gc[2, i] < old_gc[2, i]) {
#             new_gc[, i - 0:1] <- old_gc[, i - 0:1]
#         }
#     }
#     if (reset) {
#         options('.bls_gc_mem' = NULL)
#     } else {
#         options('.bls_gc_mem' = new_gc)
#     }
#     invisible(new_gc)
# }
.record_mem <- function(start = FALSE, reset = FALSE) {
    if (start) options('.bls_memory' = NULL)
    # possible alternatives:
    ps <- unlist(strsplit(system(paste0('ps -u$USER -o comm,rss,vsz,pid | grep "\\<', Sys.getpid(), '"$'), intern = TRUE), split = ' +'))
    # check process name
    stopifnot(ps[1] == 'R')
    # get RSS and vsize
    new_mem <- paste(new_num <- c(rss = as.numeric(ps[2]) / 1024, vsz = as.numeric(ps[3]) / 1024), 'MiB')
    # get old mem
    old_mem <- getOption('.bls_memory', rep('0 MiB', 2))
    old_num <- as.numeric(sub(' MiB', '', old_mem))
    for (i in 1:2) {
        if (old_num[i] > new_num[i]) {
            new_mem[i] <- old_mem[i]
        }
    }
    if (reset) {
        options('.bls_memory' = NULL)
    } else {
        options('.bls_memory' = new_mem)
    }
    invisible(c(ps[4], new_mem))
}
.start_recording <- function() {
    options(
        '.bls_record_mem' = TRUE
        )
}
.stop_recording <- function() {
    options(
        '.bls_record_mem' = NULL
        )
}
.is_recording <- function() {
    getOption('.bls_record_mem', FALSE)
}
.set_recording <- function(x) {
    options('.bls_record_mem' = isTRUE(x))
}
.record_now <- function(start = FALSE, reset = FALSE, envir = parent.frame()) {
    if (getOption('.bls_record_mem', FALSE)) {
        # return(
        #     invisible(
        #         list(
        #             gc_mem = .record_gc_mem(start, reset),
        #             object_sizes = .record_object_sizes(start, reset, envir = envir)
        #         )
        #     )
        # )
        return(.record_mem(start, reset))
    }
    invisible()
}
.gather_mem <- function(out_list) {
    mem_list <- lapply(out_list, attr, 'cpu_mem')
    mem <- rbindlist(lapply(mem_list, as.list))
    if (nrow(mem) == 0) return(invisible())
    if (names(mem)[1] == 'pid') return(mem)
    mem[, {
        num2 <- as.numeric(sub(' MiB', '', V2))
        num3 <- as.numeric(sub(' MiB', '', V3))
        out <- setNames(
                sprintf(
                    '%1.1f MiB',
                    c(
                        rss_avg = mean(num2),
                        rss_max = max(num2),
                        vsz_avg = mean(num3),
                        vsz_max = max(num3)
                    )
                ),
                c('rss_avg', 'rss_max', 'vsz_avg', 'vsz_max')
            )
        as.list(out)
    }, by = .(pid = V1)]
}

memory_usage <- function(res, show = TRUE, slurm = NULL) {
    mem <- attr(res, 'cpu_mem')
    if (is.null(mem)) return(invisible())
    num <- mem[, lapply(.SD, function(x) as.numeric(sub(' MiB', '', x)))]
    out <- num[, .(
        n_cpus = .N,
        rss_cpu = max(rss_max),
        vsz_cpu = max(vsz_max),
        rss_avg = sum(rss_avg),
        rss_max = sum(rss_max),
        vsz_avg = sum(vsz_avg),
        vsz_max = sum(vsz_max)
        )]
    fs <- function(x, base = 1024) format(structure(x * base ^ 2, class = 'object_size'),
        units = 'auto', standard = 'SI')
    has_slurm <- !is.null(slurm)
    msg <- out[, paste0(
        '~~~~~ memory usage ~~~~~\n',
        if (has_slurm) {
            paste0(
                slurm$part[, nodes], ' Nodes\n',
                n_cpus, ' CPUs (available: ', slurm$part[, cpus], ')\n'
            )
        } else {
            paste0(n_cpus, ' CPUs\n')
        },
        '\n',
        'per CPU\n',
        '=======\n',
        if (has_slurm) {
            paste0(
                'available:         ', fs(slurm$part[, total_memory / n_cpus], 1000), '\n',
                'threshold (min):   ', fs(slurm$part[, minimum_mem_given], 1000), '\n'
            )
        },
        'RSS max:           ', fs(rss_cpu), '\n',
        'VSZ max:           ', fs(vsz_cpu), '\n',
        '\n',
        'total\n',
        '=====\n',
        if (has_slurm) {
            paste0('available:         ', fs(slurm$part[, total_memory], 1000), '  (', 
                slurm$part[, nodes],' \u00D7 ', fs(slurm$part[, total_memory / nodes], 1000) ,')\n')
        },
        'RSS avg:           ', fs(rss_avg), '\n',
        'RSS max:           ', fs(rss_max), '\n',
        'VSZ avg:           ', fs(vsz_avg), '\n',
        'VSZ max:           ', fs(vsz_max), '\n',
        '~~~~~~~~~~~~~~~~~~~~~~~~\n'
        )]
    if (show) {
        cat(msg)
    }
    invisible(structure(out, msg = msg))
}
