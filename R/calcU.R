calcU <-
function(ustar,Zo,L,z,kv=0.4){
zL <- z/L
ZoL <- Zo/L
psiMz <- ifelse(zL<0,
{x <- (1-16*zL)^(1/4);-2*log((1 + x)/2) - log((1 + x^2)/2) + 2*atan(x) - pi/2},
4.8*zL
)
psiMZo <- ifelse(zL<0,
{x <- (1-16*ZoL)^(1/4);-2*log((1 + x)/2) - log((1 + x^2)/2) + 2*atan(x) - pi/2},
4.8*ZoL
)
ustar/kv*(log(z/Zo) + psiMz - psiMZo)
}
