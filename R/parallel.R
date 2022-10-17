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
