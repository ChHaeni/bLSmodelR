

print.Interval <- function(x, show_all = TRUE, n_split = 4, ...){
    x[,12] <- sapply(x[,12],function(y){
        if(!is.na(y)){
            rplc <- sub(".{0,20}$","",sub("^.{20}","",y))
        } else {
            rplc <- "" 
        }
        if(rplc!=""){
            sub(rplc,"//...//",y,fixed=TRUE)
        } else {
            y
        }
    }, USE.NAMES = FALSE)
    x[,13] <- sapply(x[,13],function(y){
        if(!is.na(y)){
            rplc <- sub(".{0,20}$","",sub("^.{20}","",y))
        } else {
            rplc <- "" 
        }
        if(rplc!=""){
            sub(rplc,"//...//",y,fixed=TRUE)
        } else {
            y
        }
    }, USE.NAMES = FALSE)
    if (!show_all) {
        if (length(x) > 13) cat(" ->", length(x) - 13, 
            "columns providing model unspecific data not shown <-\n")
        x <- x[, 1:13]
    }
    class(x) <- 'data.frame'
    if (nrow(x) > (n_split * 2)) {
        y <- rbind(x[1:n_split, ], rep(NA, ncol(x)), x[nrow(x) - ((n_split - 1):0), ])
        rownames(y) <- c(as.character(1:n_split), "", as.character(nrow(x) - ((n_split - 1):0)))
        z <- capture.output(print(y))
        lin <- grepl('^\\s+(NA|<NA>).*', z)
        z[lin] <- gsub("<NA>", " ***", z[lin])
        z[lin] <- gsub(" NA", "***", z[lin])
        z[1] <- paste0(" ", z[1])
        cat(paste0(z, "\n"))
    } else {
        print(x)
    }
}


