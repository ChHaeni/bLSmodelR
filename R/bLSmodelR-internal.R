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
.MatchWrapper <- function(Int_Ext,Cat_list,Tol_){
	# check matching
	Int_Ext[,{
		# browser()
		Cat_list[
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
			devZSens = abs(Cat_ZSens/SensorHeight - 1)/Tol_[1],
			devL = abs(Cat_L/L - 1)/Tol_[2],
			devZo = abs(Cat_Zo/Zo - 1)/Tol_[3],
			devsUu = abs(Cat_Su_Ustar/sUu - 1)/Tol_[4],
			devsVu = abs(Cat_Sv_Ustar/sVu - 1)/Tol_[5],
			devSensor_Swustar = abs(Cat_Sensor_Swustar/Sensor_Swustar - 1)/Tol_[6],
			devN0 = Cat_N0 - N0
			)
		][
			devZSens <= 1.0000001 &
			devL <= 1.0000001 &
			devZo <= 1.0000001 &
			devsUu <= 1.0000001 &
			devsVu <= 1.0000001 &
			devSensor_Swustar <= 1.0000001,sumDev:=devZSens+devL+devZo+devsUu+devsVu+devSensor_Swustar]

}
.CrossMatchWrapper <- function(Int_Ext,Cat_list,Tol_){

			# check matching heights etc.
			CheckPart1 <- Int_Ext[, .(rn, SensorHeight,Sensor_Swustar,MaxFetch,
				sWu_Upper,
				sWu_Lower,
				SensorHeight_Upper,
				SensorHeight_Lower
				)][,{
				Cat_list[
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
			Int_Ext[,.(rn,Cat.Name,SensorHeight,L,Zo,sUu,sVu,Sensor_Swustar,MaxFetch,Sensor,N0,
				sUu_Upper,
				sUu_Lower,
				sVu_Upper,
				sVu_Lower,
				Zo_Upper,
				Zo_Lower,
				L_Upper,
				L_Lower
				)][,{

				ind <- paste(SensorHeight, Sensor_Swustar, MaxFetch, sep = "/")
				Sub <- CheckPart1[.(ind), .(Zeile, ssm)]
				ind2 <- Sub[,unique(Zeile)]
				out <- Cat_list[match(ind2,Zeile)][
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
			}, by = rn][
				,":="(
				devZSens = abs(Cat_ZSens/SensorHeight - 1)/Tol_[1],
				devL = abs(Cat_L/L - 1)/Tol_[2],
				devZo = abs(Cat_Zo/Zo - 1)/Tol_[3],
				devsUu = abs(Cat_Su_Ustar/sUu - 1)/Tol_[4],
				devsVu = abs(Cat_Sv_Ustar/sVu - 1)/Tol_[5],
				devSensor_Swustar = abs(Cat_Sensor_Swustar/Sensor_Swustar - 1)/Tol_[6],
				devN0 = Cat_N0 - N0
				)
			][
				devZSens <= 1.0000001][
				devL <= 1.0000001][
				devZo <= 1.0000001][
				devsUu <= 1.0000001][
				devsVu <= 1.0000001][
				devSensor_Swustar <= 1.0000001,sumDev:=devZSens+devL+devZo+devsUu+devsVu+devSensor_Swustar]

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
