calcsigmaW <-
function(ustar,zL,bw){
phiW <- ifelse(zL<0,
(1 - 3*zL)^(1/3),
rep(1,length(zL))
)
return(bw*ustar*phiW)
}
