round125 <- function (x, rx = c(1, 2, 5)){
    rx2 <- c(rx, rx[1] * 10)
    rxm <- c(0, (rx + rx2[-1])/2)
    lo10 <- 10^floor(log10(abs(x)))
    m <- 1/lo10
    y <- x * m * sign(x)
    ind <- sapply(y, function(y) rev(which(rxm <= y))[1])
    out <- rx2[ind]/m * sign(x)
    return(out)
}


