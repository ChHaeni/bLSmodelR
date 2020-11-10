calcbw <- function(sigmaWustar,zL){
	phiW <- ifelse(zL<0,
		(1 - 3*zL)^(1/3),
		rep(1,length(zL))
	)
	return(sigmaWustar/phiW)
}
