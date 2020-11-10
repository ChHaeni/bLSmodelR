genTolerances <- function(Tol.Z=2.5,Tol.L=10,Tol.Zo=5,Tol.sUu=10,Tol.sVu=10,Tol.sWu=5,Tol.Zero=FALSE){
	if(Tol.Zero){
		Out <- data.frame("Tolerance (in %)"=rep(0,9),stringsAsFactors=FALSE,check.names=FALSE)
	} else {
		Out <- data.frame("Tolerance (in %)"=c(Tol.Z,Tol.L,Tol.Zo,Tol.sUu,Tol.sVu,Tol.sWu,0,0,0),stringsAsFactors=FALSE,check.names=FALSE)
	}
	if(any(is.na(Out)))stop("NA values are not allowed...\n")
	Out <- round(Out/5,1)*5
	rownames(Out) <- c("Sensor Height","L","Zo","SigmaU/Ustar","SigmaV/Ustar","SigmaW/Ustar","kv","A","alpha")
	return(structure(Out,class=c("Tolerances","data.frame")))
}
