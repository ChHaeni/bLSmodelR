calcepsilon <- function(ustar,L,bw,z,kv=0.4){
	zL <- z/L
	phiEpsilon <- ifelse(zL<0,
		(bw^4*(1 - 3*zL)^(4/3) + 1)/((bw^4 + 1)*(1 - 3*zL)^(1/3)*(1 - 6*zL)^(1/4)),
		1 + 5*zL
	)

	return(ustar^3*phiEpsilon/(kv*z))
}
