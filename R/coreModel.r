coreModel <- function(u, v, w, zSens, ustar, L, Zo, bw, sigmaUustar, sigmaVustar, kv, C0, alpha, MaxFetch){
    if (add <- any(inherits(u, 'TDcat'))) {
		name = substitute(u)
	    Ctlg <- u
		Head <- unlist(strsplit(attr(u,"header"),"\n"))[-1]
		Whead <- matrix(as.numeric(gsub(".*[=] ", "", Head)), nrow = 1)
		zSens <- Whead[,2]
		ustar <- Whead[,3]
		L <- Whead[,4]
		Zo <- Whead[,5]
		sigmaUustar <- Whead[,6]
		sigmaVustar <- Whead[,7]
		bw <- Whead[,8]
		C0 <- Whead[,9]
		kv <- Whead[,10]
		alpha <- Whead[,12]
		MaxFetch <- Whead[,13]

		v <- attr(u,"uvw0")[,2]
		w <- attr(u,"uvw0")[,3]
		u <- attr(u,"uvw0")[,1]
	}
	Linv <- 1/L
    if (Linv < 0) {
        out <- csFi(u, v, w, zSens, ustar, Linv, Zo, bw, sigmaUustar, 
            sigmaVustar, kv, C0, alpha, MaxFetch)
    } else {
        out <- csFs(u, v, w, zSens, ustar, Linv, Zo, bw, sigmaUustar, 
            sigmaVustar, kv, C0, alpha, MaxFetch)
    }
	if(add){
		attCat <- attributes(Ctlg)
	    Ctlg <- rbindlist(list(Ctlg, out))
	    for (ac in (names(attCat) %w/o% c("names", "row.names", ".internal.selfref"))) setattr(Ctlg, 
	        ac, attCat[[ac]])
	    assign(as.character(name), Ctlg, parent.frame(), inherits = TRUE)
	    return(invisible(Ctlg))
	} else {
	    return(out)
	}
}


