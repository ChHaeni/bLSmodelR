# functions to fix bug in parallel::clusterApplyLB
.clusterApplyLB <- function (cl = NULL, x, fun, ...) {
    argfun <- function(i) c(list(x[[i]]), list(...))
    .dynamicClusterApply(cl, fun, length(x), argfun)
}
.dynamicClusterApply <- function (cl = NULL, fun, n, argfun) {
    cl <- parallel:::defaultCluster(cl)
    p <- length(cl)
    if (n > 0L && p) {
        submit <- function(node, job) parallel:::sendCall(cl[[node]], fun, 
            argfun(job), tag = job)
        for (i in 1:min(n, p)) submit(i, i)
        val <- vector("list", n)
        for (i in 1:n) {
            d <- .recvOneResult(cl)
            j <- i + min(n, p)
            if (j <= n) 
                submit(d$node, j)
            val[d$tag] <- list(d$value)
        }
        parallel:::checkForRemoteErrors(val)
    }
}
.recvOneResult <- function (cl) {
    if (parallel:::.snowTimingData$running()) {
        start <- proc.time()[3]
        v <- .recvOneData_SOCK(cl)
        end <- proc.time()[3]
        parallel:::.snowTimingData$enterRecv(v$node, start, end, v$value$time[3])
    }
    else v <- .recvOneData_SOCK(cl)
    list(value = v$value$value, node = v$node, tag = v$value$tag)
}
.recvOneData_SOCK <- function (cl) {
    socklist <- lapply(cl, function(x) x$con)
    repeat {
        ready <- socketSelect(socklist)
        if (any(ready)) 
            break
    }
    n <- which(ready)[1L]
    list(node = n, value = unserialize(socklist[[n]]))
}

### fix memory limit
.makePSOCKcluster <- function (names, memory_limit = NA, ...) {
    options <- parallel:::addClusterOptions(parallel:::defaultClusterOptions, list(...))
    manual <- parallel:::getClusterOption("manual", options)
    homogeneous <- parallel:::getClusterOption("homogeneous", options)
    setup_strategy <- match.arg(parallel:::getClusterOption("setup_strategy", 
        options), c("sequential", "parallel"))
    setup_timeout <- parallel:::getClusterOption("setup_timeout", options)
    local <- is.numeric(names) || (is.character(names) && identical(names, 
        rep("localhost", length(names))))
    if (is.numeric(names)) {
        names <- as.integer(names[1L])
        if (is.na(names) || names < 1L) 
            stop("numeric 'names' must be >= 1")
        names <- rep("localhost", names)
    }
    parallel:::.check_ncores(length(names))
    cl <- vector("list", length(names))
    if (!manual && homogeneous && local && setup_strategy == 
        "parallel") {
        port <- parallel:::getClusterOption("port", options)
        timeout <- parallel:::getClusterOption("timeout", options)
        useXDR <- parallel:::getClusterOption("useXDR", options)
        cmd <- parallel:::workerCommand("localhost", options, setup_strategy = "parallel")
        socket <- serverSocket(port = port)
        on.exit(close(socket), add = TRUE)
        if (.Platform$OS.type == "windows") {
            for (i in seq_along(cl)) system(cmd, wait = FALSE, 
                input = "")
        }
        else {
            #~~~~ memory limit
            if (!missing(memory_limit))
                cmd <- paste0('export R_MAX_VSIZE=', memory_limit, ' && ', cmd)
            #~~~~ memory limit
            cmd <- paste(rep(cmd, length(cl)), collapse = " & ")
            system(cmd, wait = FALSE)
        }
        cls <- if (useXDR) 
            "SOCKnode"
        else "SOCK0node"
        ready <- 0
        pending <- list()
        on.exit(lapply(pending, function(x) close(x$con)), add = TRUE)
        t0 <- Sys.time()
        while (ready < length(cl)) {
            cons <- lapply(pending, function(x) x$con)
            if (difftime(Sys.time(), t0, units = "secs") > setup_timeout + 
                5) {
                failed <- length(cl) - ready
                msg <- sprintf(ngettext(failed, "Cluster setup failed. %d worker of %d failed to connect.", 
                  "Cluster setup failed. %d of %d workers failed to connect."), 
                  failed, length(cl))
                stop(msg)
            }
            a <- socketSelect(append(list(socket), cons), FALSE, 
                timeout = setup_timeout)
            canAccept <- a[1]
            canReceive <- seq_along(pending)[a[-1]]
            if (canAccept) {
                con <- socketAccept(socket = socket, blocking = TRUE, 
                  open = "a+b", timeout = timeout)
                scon <- structure(list(con = con, host = "localhost", 
                  rank = ready), class = cls)
                tryCatch({
                    parallel:::sendCall(scon, eval, list(quote(Sys.getpid())))
                }, error = identity)
                pending <- append(pending, list(scon))
            }
            for (scon in pending[canReceive]) {
                pid <- tryCatch({
                    parallel:::recvResult(scon)
                }, error = identity)
                if (is.integer(pid)) {
                  ready <- ready + 1
                  cl[[ready]] <- scon
                }
                else close(scon$con)
            }
            if (length(canReceive) > 0) 
                pending <- pending[-canReceive]
        }
    }
    else {
        for (i in seq_along(cl)) cl[[i]] <- parallel:::newPSOCKnode(names[[i]], 
            options = options, rank = i)
    }
    class(cl) <- c("SOCKcluster", "cluster")
    cl
}


