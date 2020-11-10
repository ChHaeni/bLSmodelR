genModel <- function(kv=0.4,A=0.5,alpha=0.02,wTDcutoff=1E-4,TDwrite=TRUE,overwriteTD=TRUE,TDread=TRUE,TDonly=FALSE,ncores=1){
	Out <- list(kv,A,alpha,wTDcutoff,TDwrite,overwriteTD,TDread,TDonly,ncores)
	names(Out) <- nams <- c("kv","A","alpha","wTDcutoff","TDwrite","overwriteTD","TDread","TDonly","ncores")
	if(any(is.na(Out)))stop("NA values are not allowed...\n")
	if(any(wn <- !sapply(Out[c(1:4,9)],is.numeric)))stop("Argument ",paste0(nams[c(1:4,9)][wn],collapse=", ")," must be numeric!\n")
	if(any(wn <- !sapply(Out[5:8],is.logical)))stop("Argument ",paste0(nams[5:8][wn],collapse=", ")," must be logical!\n")
	Out[[1]] <- round(Out[[1]],2)
	Out[[2]] <- round(Out[[2]],2)
	Out[[3]] <- round(Out[[3]],3)
	if(Out[[3]]<1E-3){
		warning("alpha value is too small. value has been reset to 0.001!")
		Out[[3]] <- 1E-3
	}
	if(Out[[4]]<1E-4){
		warning("wTDcutoff value is too small. value has been reset to 1E-4!")
		Out[[4]] <- 1E-4
	}
	Out[[9]] <- as.integer(Out[[9]])
	return(structure(Out,class=c("Model","list")))
}
