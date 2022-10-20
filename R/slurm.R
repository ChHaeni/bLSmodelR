# runbLSlurm main function
runbLSlurm <- function(input_list, cat_path, ..., 
    memory_limit = NULL, record_mem = FALSE, wait = TRUE) {
    
    # print usage without arguments
    if ((missing(input_list) || missing(cat_path) || missing(...))) {
        # print usage here (especially what to pass as ...)
        cat('"runbLSlurm" arguments:
        \r\tinput_list - list with model input (ModelInput in runbLS)
        \r\tcat_path   - path to catalog directory (Cat.Path in runbLS)
        \r\t...        - slurm options which are transfered to the sbatch file
        \r\t             mandatory is "mem". 
        \r\t             "partition" and "nodes" can be provided. 
        \r\t             If they are not provided, the best partition/node/cpu combination will be chosen.
        \r\t             "cpus", resp. "cpus-per-task" is derived from the number of rows that need to be calculated.
        \r\t             These arguments can also be provided as list().
        \r\twait       - if TRUE (default), the function will wait for the slurm job to finish
        \r\t             and collect the slurm job results before returning it.
        \r\t             If FALSE, the job-id and the job-directory will be returned such that the user
        \r\t             can collect the results via the function "collect_results".\n')
        return(invisible(NULL))
    }

    # memory recording?
    options('.bls_record_mem' = isTRUE(record_mem))
    on.exit(options('.bls_record_mem' = NULL), add = TRUE)

    # remove NA values in Interval data.frame
    isna <- as.logical(rowSums(is.na(input_list$Interval[, 1:13])))
    input_list$Interval <- input_list$Interval[!isna, ]
    if (any(isna)) cat('Removed', sum(isna), 'rows due to NA values\n')

    # extract slurm options and prepare job directory
    slurm <- prep_slurm(..., ntasks = NROW(input_list$Interval))

    # split Intervals and save to rds files
    il <- split_int(input_list$Interval, slurm$part)
    for (i in seq_along(il)) {
        saveRDS(il[[i]], file.path(slurm$tmp_dir, paste0('int', i, '.rds')))
    }

    # remove Interval and save model input list
    input_list$Interval <- NULL
    saveRDS(input_list, file.path(slurm$tmp_dir, 'input_list.rds'))

    # create script with argument
    rscript_file <- write_runbLS_script(slurm$tmp_dir, cat_path, slurm$part[, cpus_per_task], mem_lim = memory_limit)

    # create sbatch file, run slurm job & return result
    run_sbatch(slurm = slurm, rscript = rscript_file, wait = wait)

}

# create and run sbatch file
run_sbatch <- function(slurm, rscript, wait) {
    # dots without partition, nodes, cpus_per_task
    dots <- get_sopt(slurm$dots, 
        'p', 'partition', 'N', 'nodes', 'cpus.*', 'J', 'job-name', 
        remove = TRUE) 
    # create sbatch file
    sbatch_file <- write_sbatch(
        slurm$tmp_dir, 
        rscript,
        'job-name' = slurm$job_name,
        partition = slurm$part[, sub('[*]', '', Part)],
        nodes = slurm$part[, nodes],
        'cpus-per-task' = slurm$part[, cpus_per_task],
        dots
    )
    # run sbatch file
    re <- system(paste('sbatch', sbatch_file), intern = TRUE)
    cat(re, '\n')
    # capture job id
    job_id <- sub('[^0-9]*([0-9]*)[^0-9]*', '\\1', re)
    # if else: either wait for job to finish or save job id and collect later
    if (wait) {
        # get current time
        current_time <- Sys.time()
        # wait for job to finish
        # how can we know that jobs finished?
        # squeue -j $job_id
        sq <- system(paste0('squeue -j ', job_id), intern = TRUE)
        # be verbose
        cat(sq[1], '\n')
        cat(sq[2])
        # initiate seconds to wait
        sleep_secs <- 1
        # loop until job is done
        while (length(sq) == 2) {
            # wait some time
            Sys.sleep(sleep_secs)
            # get job state
            sq <- system(paste0('squeue -j ', job_id), intern = TRUE)
            if (length(sq) == 2) {
                cat('\r', sq[2])
                # TODO: grep for job state and start time recording once job started
            }
            # increase waiting time by 1 second
            sleep_secs <- sleep_secs + 1
        }
        cat('\njob finished.\n')
        # job summary
        cat('\n')
        seff(job_id)
        cat('\n')
        # collect and return results
        res <- collect_results(slurm$tmp_dir)
        # duration of job?
        dur <- Sys.time() - current_time
        cat('Time since sending job: ', round(dur, 2), attr(dur, 'units'), '\n')
        # return value
        res
    } else {
        # be verbose
        cat('Not waiting for job to finish & returning job info.
        \rUse "collect_results(x)" where argument "x" is either
        \r\tthe returned object or
        \r\tthe path to the job directory
        \r')
        # return job id and slurm$tmp_dir
        list(
            'job-dir' = slurm$tmp_dir,
            'job-id' = job_id,
            'job-name' = slurm$job_name,
            'partition' = slurm$part[, Part],
            'nodes' = slurm$part[, nodes],
            'nodelist' = slurm$part[, node_names],
            'cpus-per-task' = slurm$part[, cpus_per_task]
        )
    }
}

# write sbatch file function
write_sbatch <- function(tmpdir, rscript, ...) {
    # have a look at sample slurm scripts: https://help.rc.ufl.edu/doc/Sample_SLURM_Scripts
    # and here!: https://support.ceci-hpc.be/doc/_contents/QuickStart/SubmittingJobs/SlurmTutorial.html
    # capture sbatch arguments
    arg_values <- list(...)
    # check if ... has list
    if (any(sapply(arg_values, is.list))) {
        arg_values <- unlist(arg_values, recursive = FALSE)
    }
    # check for argument names
    arg_names <- names(arg_values)
    # error if unnamed argument(s) exist
    if (is.null(arg_names) || any(arg_names %in% '')) {
        stop('Unnamed SBATCH arguments')
    }
    # get temporary file name
    tmp <- tempfile(pattern = 'sbatch', tmpdir = tmpdir)
    # find Rscript path
    rscript_path <- file.path(R.home(), 'bin/Rscript')
    # write to file
    writeLines(c(
        '#!/bin/bash',
        # change working directory
        paste0('#SBATCH --chdir=', tmpdir),
        # loop over arguments
        sapply(arg_names, function(a) {
            # remove dash if nchar is 2
            if (nchar(a) == 2) a <- sub('[-]', '', a)
            # check argument
            arg <- switch(as.character(nchar(a))
                # add '-' if one char
                , '1' = paste0('-', a)
                # prepend '--'
                , sub('^[-]{0,2}', '--', a)
            )
            # concatenate string
            paste0('#SBATCH ', arg, '=', arg_values[[a]])
        }),
        paste0('for file in ', tmpdir, '/int*.rds\ndo'),
        paste0('srun -n1 -N1 ', rscript_path, ' --vanilla ', rscript, ' $file &\ndone\nwait')
        ),
        tmp
    )
    # return sbatch file path
    tmp
}

# prepare slurm arguments and directory
prep_slurm <- function(..., ntasks = 1) {
    # get dot arguments
    dots <- list(...)
    # check if list of options has been provided
    if (any(sapply(dots, is.list))) {
        dots <- unlist(dots, recursive = FALSE)
    }
    # remove -n --ntasks and give warning
    dots <- clean_ntasks(dots)
    # get current time
    current_time <- Sys.time()
    # jobname (-J --job-name)
    job_name <- get_jobname(dots, current_time)
    # memory usage (--mem --mem-per-cpu)
    #    an option could be, to distribute mem and tasks
    #    to different nodes with different cpus
    #    for now, finding the most cpus by using equal mem/cpu
    #    distribution should be sufficient
    mem <- get_mem(dots)
    # find partition and nodes
    #    partition name (-p --partition)
    #    number of nodes (-N --nodes)
    #    cpus per task
    part <- find_partition(mem, ntasks, dots)
    # create temporary directory
    #    check if -D --chdir exists
    #    otherwise default to $HOME/.slurm/$jobname
    tmp_dir <- get_sopt(dots, 'D', 'chdir', alternative = {
        # get directory from options
        job_parent_dir <- getOption(
            'bls.slurm.jobdir',
            default = file.path(Sys.getenv('HOME'), '.slurm')
        )
        # set path to $HOME/.slurm
        file.path(job_parent_dir, job_name)
    })
    # create tmp_dir if directory doesn't exists
    if (!dir.exists(tmp_dir)) {
        dir.create(tmp_dir, recursive = TRUE)
    } else {
        cat(paste0('Directory "', tmp_dir, '" already exists with content:\n'))
        system(paste('ls -la', tmp_dir))
        cat('Make sure that you have removed any unwanted files in it!\n')
        ans <- 'ask'
        while (!(ans %in% c('y', 'Y', 'yes', 'N', 'n', 'no', ''))) {
            ans <- readline('Do you want to proceed? [Y/n]: ')
        }
        if (any(ans %in% c('N', 'n', 'no'))) return(invisible(NULL))
    }
    # return list with dots & part
    list(
        dots = dots,
        part = part,
        tmp_dir = tmp_dir,
        job_name = job_name
        )
}

# split Interval into chunks
split_int <- function(int, p){
    # get rows
    nr <- NROW(int)
    # get nodes
    nodes <- p[, nodes]
    # minimum of tasks per node
    mt <- floor(nr / nodes)
    # residual tasks
    rt <- nr %% nodes
    # number of tasks per node
    nt <- c(rep(mt + 1, rt), rep(mt, nodes - rt))
    # upper
    up <- cumsum(nt)
    # lower
    lo <- c(1, up[-length(up)] + 1)
    # index
    ind <- mapply(':', lo, up, SIMPLIFY = FALSE)
    # split
    lapply(ind, function(i) int[i, ])
}

# write R scripts
# runbLS
write_runbLS_script <- function(tmpdir, cpath, ncores, mem_lim = NULL) {
    # get tmpfile name
    tmp <- tempfile(pattern = 'Rscript', tmpdir = tmpdir, fileext = '.R')
    # write R script to tmp file
    writeLines(
        c(
            'library(bLSmodelR)',
            # memory recording?
            paste0('options(".bls_record_mem" = ', getOption('.bls_record_mem', FALSE)),
            # format of file: int%i.rds
            'ifile <- commandArgs(TRUE)',
            # read intervals
            'int <- readRDS(ifile)',
            # read inlist
            paste0('inlist <- readRDS(file.path("', tmpdir, '", "input_list.rds"))'),
            # add int to inlist
            'inlist$Interval <- int',
            # run model
            if (is.null(mem_lim)) {
                paste0('res <- runbLS(inlist, "', cpath, '", ncores = ', ncores, ')')
            } else {
                paste0('res <- runbLS(inlist, "', cpath, '", ncores = ', ncores, ', memory_limit = "', mem_lim, '")')
            },
            # save result; get index from int%i.rds
            'saveRDS(res, sub("/int([0-9]{1,2}[.]rds)", "/res\\\\1", ifile))'
        ), 
        tmp
    )
    # make file executable
    # Sys.chmod could be a better option?
    system(paste('chmod +x', tmp))
    # return tmpfile name
    tmp
}
# depostion
write_deposition_script <- function(tmpdir, ncores, mem_lim = NULL) {
    # get tmpfile name
    tmp <- tempfile(pattern = 'Rscript', tmpdir = tmpdir, fileext = '.R')
    # write R script to tmp file
    writeLines(
        c(
            'library(bLSmodelR)',
            # memory recording?
            paste0('options(".bls_record_mem" = ', getOption('.bls_record_mem', FALSE)),
            # format of file: int%i.rds
            'ifile <- commandArgs(TRUE)',
            # read intervals
            'int <- readRDS(ifile)',
            # read bls result
            'bls_result <- int',
            # read arguments and attributes
            paste0('dep_args <- readRDS(file.path("', tmpdir, '", "dep_args.rds"))'),
            # add attributes to result
            'attr(bls_result, "ModelInput") <- dep_args[["ModelInput"]]',
            'attr(bls_result, "Catalogs") <- dep_args[["Catalogs"]]',
            'attr(bls_result, "Cat.Path") <- dep_args[["Cat.Path"]]',
            # run deposition
            if (is.null(mem_lim)) {
                paste0('res <- do.call(deposition, c(list(x = bls_result, ncores = ', 
                    ncores, '), dep_args[c("vDep", "vDepSpatial")]))')
            } else {
                paste0('res <- do.call(deposition, c(list(x = bls_result, ncores = ', 
                    ncores, ', memory_limit = "', mem_lim, '"), dep_args[c("vDep", "vDepSpatial")]))')
            },
            # save result; get index from int%i.rds
            'saveRDS(res, sub("/int([0-9]{1,2}[.]rds)", "/res\\\\1", ifile))'
        ), 
        tmp
    )
    # make file executable
    # Sys.chmod could be a better option?
    system(paste('chmod +x', tmp))
    # return tmpfile name
    tmp
}

depoSlurm <- function(x, vDep, ..., rn = NULL, Sensor = NULL, Source = NULL, vDepSpatial = NULL,
    memory_limit = NULL, record_mem = FALSE, wait = TRUE) {

    # memory recording?
    options('.bls_record_mem' = isTRUE(record_mem))
    on.exit(options('.bls_record_mem' = NULL), add = TRUE)

    # convert old versions 
    sx <- as.character(substitute(x))
    x <- copy(x)
    setDT(x)
    switchNames(x)
    if(is.null(attr(x, "Version"))){
        warning(paste0("Object '", sx[min(length(sx), 2)], "' has not yet been converted to version 4.2+"))
        convert(x)
    }
    # check rn argument
	if(is.null(rn)){
		Selrn <- x[,unique(rn)]
	} else {
		Selrn <- rn
	}
    # check Sensor argument
	if(is.null(Sensor)){
		SelSensor <- x[,unique(Sensor)]
	} else {
		SelSensor <- Sensor
	}
    # check Source argument
	if(is.null(Source)){
		SelSource <- x[,unique(Source)]
	} else {
		SelSource <- Source
	}
    # subset
	Run <- x[rn %in% Selrn & Sensor %chin% SelSensor & Source %chin% SelSource,]

    # throw error if vDep/vDepSpatial length > 1
    if (length(vDep) != 1) {
        stop('Argument "vDep" must be either a single number or a column name')
    }
    if (!is.null(vDepSpatial) && (!is.list(vDepSpatial) || 
        length(vDepSpatial) != 2 || any(!(names(vDepSpatial[[1]]) %in% unique(vDepSpatial[[2]][, 1]))))) {
        stop('Argument "vDepSpatial" must be a list with two elements:
            \r\t first element: named list with corresponding "extra" vDep values (or column names)
            \r\t second element: data.frame of class "Source" with "extra" areas
            \r\t               -> any specified "extra" name in vDepSpatial[[1]] must have a
            \r\t                  corresponding source element in vDepSpatial[[2]]!'
        )
    }

    # extract slurm options and prepare job directory
    slurm <- prep_slurm(..., ntasks = Run[, .N])

    # split Intervals and save to rds files
    il <- split_int(Run, slurm$part)
    for (i in seq_along(il)) {
        # strip attributes
        attributes(il) <- attributes(il)[c('row.names', 'class', '.internal.selfref', 'names')]
        # save to file
        saveRDS(il[[i]], file.path(slurm$tmp_dir, paste0('int', i, '.rds')))
    }

    # save function arguments and bls result attributes to rds file
    saveRDS(
        list(
            # arguments
            vDep = vDep, 
            vDepSpatial = vDepSpatial,
            # attributes
            ModelInput = attr(x, "ModelInput"),
            Catalogs = attr(x, "Catalogs"),
            Cat.Path = attr(x, "CatPath")
            ), 
        file.path(slurm$tmp_dir, "dep_args.rds")
    )
    
    # create script with argument
    rscript_file <- write_deposition_script(slurm$tmp_dir, slurm$part[, cpus_per_task], mem_lim = memory_limit)

    # create sbatch file, run slurm job & return result
    run_sbatch(slurm = slurm, rscript = rscript_file, wait = wait)

}

# get slurm option
find_sopt <- function(x, ...) {
    if (length(x)) {
        # loop over arguments
        out <- sapply(list(...), function(x) {
            # remove trailing dashes
            if (nchar(x) == 2) {
                x <- sub('[-]', '', x)
            } else {
                x <- sub('^[-]{2}', '', x)
            }
            # add optional dashes
            switch(as.character(nchar(x))
                , '1' = paste0('^[-]?', x, '$')
                , paste0('^([-]{2})?', x, '$')
            )
        })
        # paste with OR |
        rexpr <- paste0('(', paste(out, collapse = '|'), ')')
        grepl(rexpr, names(x))
    } else {
        NULL
    }
}

# get slurm option value
get_sopt <- function(x, ..., alternative = x, remove = FALSE) {
    ind <- find_sopt(x, ...)
    if (any(ind)) {
        if (remove) {
            x[!ind]
        } else {
            if (sum(ind) > 1) {
                x[ind]
            } else {
                x[[which(ind)]]
            }
        }
    } else {
        alternative
    }
}

# get job name function
get_jobname <- function(x, ctime) {
    jn <- find_sopt(x, 'J', 'job-name')
    if (any(jn)) {
        x[[which(jn)]]
    } else {
        format(ctime, 'job_%y%m%d_%H%M%S')
    }
}

# get memory option
get_mem <- function(x) {
    mn <- find_sopt(x, 'mem')
    if (any(mn)) {
        x[[which(mn)]]
    } else {
        stop('Please provide the amount of memory to allocate for each node (--mem=..?)')
    }
}

clean_ntasks <- function(x) {
    tn <- find_sopt(x, 'n', 'ntasks')
    if (any(tn)) {
        warning('argument "', names(x)[tn], '" will be ignored since the number of tasks is given by the number of "Interval" rows')
        x[[-which(tn)]]
    } else {
        x
    }
}

# find partition
# TODO: switch to --mem-per-cpu !!! and export seff
find_partition <- function(memory, ...) {
    # ni call
    ni_call <- 'ni'
    # capture dots
    dts <- list(...)
    # check if list of options has been provided
    if (any(sapply(dts, is.list))) {
        dts <- unlist(dts, recursive = FALSE)
    }
    # print table without arguments
    if (missing(memory) && length(dts) == 0) {
        return(system('ni'))
    } else if (is.na(suppressWarnings(as.numeric(sub('[A-Z]$', '', memory))))) {
        pattern <- paste0('.*', memory, '.*|$')
        return(system(paste0('ni | grep --color -E \'', pattern, '\'')))
    }
    # check ntasks
    ntasks <- get_sopt(dts, 'n', 'ntasks', alternative = Inf)
    # correct memory
    if (is.character(memory)) {
        mem <- suppressWarnings(as.numeric(sub('(T|G|M|K)$', '', memory)))
        memory <- mem * switch(sub('.*(T|G|M|K)$', '\\1', memory)
            , 'T' = 1e6
            , 'G' = 1e3
            , 'M' = 1
            , 'K' = 1e-3
            , stop('Memory unit not recognized')
        )
    }
    # check for partition in options
    part <- get_sopt(dts, 'p', 'partition', alternative = NULL)
    # check for number of nodes in partition
    Nodes <- get_sopt(dts, 'N', 'nodes', alternative = NULL)
    # check if numeric
    if (!is.null(Nodes)) {
        suppressWarnings(Nodes <- as.integer(Nodes))
        if (is.na(Nodes)) {
            stop('slurm option -n/--nodes must be an integer')
        }
    }
    # read partition table
    ni <- data.table::fread(cmd = ni_call)
    # any partitions to exclude?
    exclude_partitions <- getOption('bls.slurm.exclude.partition', '')
    # check if it is a vector
    if (!is.character(exclude_partitions) || length(exclude_partitions) == 0) {
        stop('option "bls.slurm.exclude.partition" should be a character vector!')
    }
    ni <- ni[!(Part %in% exclude_partitions)]
    # subset part if not null
    if (!is.null(part)) {
        ni <- ni[Part %in% part]
        if (nrow(ni) == 0) {
            stop({find_partition(); paste0('partition ', part, ' is not available')})
        }
    }
    # check memory
    ni_mem <- ni[(State %in% c('idle', 'mix')) & (MFree * 1e3) >= memory]
    if (nrow(ni_mem) == 0) {
        if (is.null(part)) {
            stop({find_partition(); paste0('No partition with enough memory available.')})
        } else {
            stop({find_partition(); paste0('Partition ', part, ' doesn\'t have nodes with enough memory available.')})
        }
    }
    # summarize
    ni_sum <- ni_mem[, {
        Cav <- unique(CIdle)
        rbindlist(lapply(Cav, function(cav) {
            # print node names if colored below
            ind <- CIdle >= cav
            Nodes <- 1:sum(ind)
            .(
                node_names = list(Node[ind]),
                nodes = Nodes,
                cpus = Nodes * cav,
                cpus_per_task = cav
            )
        }))
    }, by = Part][nodes <= ntasks]
    # check partition
    if (is.null(part)) {
        # exclude alloc
        out <- ni_sum[order(cpus, cpus_per_task, decreasing = TRUE)[1], ]
    } else {
        # select part
        out <- ni_sum[Part %in% part]
        # check if partition is available
        if (nrow(out) == 0) {
            stop({find_partition(part); paste0('Specified partition "', part, '" has not enough resources available.')})
        }
        # nodes specified or not?
        if (is.null(Nodes)) {
            out <- out[which.max(cpus), ]
        } else {
            # check if nodes number is equal
            out <- out[, {
                if (any(ind <- nodes == Nodes)) {
                    .SD[ind][which.max(cpus_per_task),]
                } else {
                    # else error
                    stop({find_partition(part); paste0('Specified partition "', part, '" has not specified number of nodes available.')})
                }
            }]
        }
    }
    # print system call with selected nodes highlighted
    patterns <- out[, paste(c(paste0(unlist(node_names), '.*'), '$'), collapse = '|')]
    system(paste0(ni_call, ' | grep --color -E \'', patterns, '\''))
    # return selected as data.table
    out
}

collect_results <- function(job_dir, check.res = TRUE) {
    # job_dir as list with jobid
    if (is.list(job_dir)) job_dir <- job_dir[['job-dir']]
    # be verbose about collecting from path and job
    cat('Collecting results from:', job_dir, '\n')
    # get res.*rds file paths
    res_files <- dir(job_dir, pattern = 'res.*[.]rds', full.names = TRUE)
    if (check.res) {
        # get int.*rds file paths
        int_files <- dir(job_dir, pattern = 'int.*[.]rds', full.names = TRUE)
        if (length(res_files) != length(int_files)) {
            stop('Some jobs didn\'t finish!\n',
                'If you still want to collect the results, run:\n',
                'collect_results("', job_dir, '", check.res = FALSE)\n'
                )
        }
    }
    # read in
    res_list <- lapply(res_files, readRDS)
    # check them
    if (any(lengths(res_list) == 0)) {
        # which empty
        ind <- which(lengths(res_list) == 0)
        # shorten list and give warning
        warning(length(ind), " job result(s) were removed because their returned values were NULL!")
        # remove empty
        res_list <- res_list[-ind]
    } else if(any(sapply(res_list, nrow) == 0)) {
        # which empty
        ind <- which(sapply(res_list, nrow) == 0)
        # shorten list and give warning
        warning(length(ind), " job result(s) were removed because their returned data.table were empty!")
        # remove empty
        res_list <- res_list[-ind]
    }
    # join them to one
    res <- do.call(join, res_list)
    # get original rn values
    rn_values <- res[, sub('^[0-9]+_', '', rn)]
    # get new rn values
    names(rn_values) <- res[, rn]
    # restore rn values
    res[, rn := rn_values]
    if (inherits(res, 'deposition')) {
        attr(res, 'vDep') <- list(
            vDep = unlist(attr(res, 'vDep')[['vDep']]),
            vDepSpatial = unlist(attr(res, 'vDep')[['vDepSpatial']])
            )
    } else {
        # also in attributes
        attr(res, 'CalcSteps')[, rn := rn_values[rn]]
        attr(res, 'Catalogs')[, rn := rn_values[rn]]
        # set keys again
        setkey(attr(res, 'CalcSteps'), rn, Sensor)
        setkey(attr(res, 'Catalogs'), rn, Sensor, PointSensor)
    }
    # return
    res
}

seff <- function(job_id, colorize = TRUE) {
    seff_out <- system(paste('seff', job_id), intern = TRUE)
    if (colorize && length(attr(seff_out, 'status')) == 0) {
        # check state COMPLETED(?)
        state <- sub('State: ([a-zA-Z0-9_.-]+) .*$', '\\1', seff_out[4])
        if (state == 'COMPLETED') {
            state <- '\033[38;5;82m~~~ SLURM JOB COMPLETED ~~~\033[0m'
            fin <- '\033[38;5;82m~~~~~~~~~~~~~~~~~\033[0m'
        } else {
            state <- paste0('\033[38;5;196m~~~ SLURM JOB ', state, ' ~~~\033[0m')
            fin <- '\033[38;5;196m~~~~~~~~~~~~~~~~~\033[0m'
        }
        cat(c(state, seff_out, fin), sep = '\n')
    } else {
        cat(seff_out, sep = '\n')
    }
    invisible(seff_out)
}



