print.Sources <- function(x,Nrows=30,...){
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
		if(nrow(x)>Nrows){
			hx <- x[seq(3),,drop=FALSE]
			tx <- x[seq.int(to = nrow(x), length.out = 3),,drop=FALSE]
			px <- data.frame(apply(rbind(hx,rep(NA,4),tx),2,as.character),check.names=FALSE)
			rownames(px) <- c(rownames(x)[1:3],"***",rownames(x)[seq.int(to = nrow(x), length.out = 3)])
			print.data.frame(px,na.print="***",quote=FALSE,...)
		} else {
			print.data.frame(x,...)
		}
	} else {
		print.data.frame(x,...)
	}
}
