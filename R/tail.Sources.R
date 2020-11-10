tail.Sources <- function(x,...){
	if(length(x)==4){
		x <- procSources(x)
		uS <- unique(x[,1])
		nS <- length(uS)
		cat("******\nSources data.frame:\n******\n")
		cat(" Number of Sources:",nS,"\n")
		SourceList <- attr(x,"SourceList")
		for(i in seq(nS)){
			p <- SourceList[[uS[i]]]
			np <- length(unique(p[,3]))
			cat(" - ",uS[i],":",round(sum(by(p,p[,3],getArea)),2),"(m2)",if(np>1){paste0("(",length(unique(p[,3]))," polygons)")},"\n")
		}
		cat("******\n\n")
	}
	class(x) <- "data.frame"
	tail(x,...)
}
