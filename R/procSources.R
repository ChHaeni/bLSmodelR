procSources <- function(Sou){
	Snames <- unique(Sou[,"SourceArea Name"])
	Slist <- lapply(Snames,function(x,y)y[y[,1]==x,2:4],Sou[,c("SourceArea Name","x-Coord (m)","y-Coord (m)","Polygon ID")])
	names(Slist) <- Snames
	SAreas <- sapply(Slist,function(x)sum(by(x,x[,3],getArea)))

	# checkSources (re-write someday!):
	for(i in 1:length(Slist)){
		Neu <- Slist[[i]][numeric(0),]
		polyIndex <- Slist[[i]][,3]
		polyU <- unique(polyIndex)
		for(j in 1:length(polyU)){
			x <- Slist[[i]][polyIndex==polyU[j],]
			if(!identical(x[1,],x[nrow(x),])){
					x <- rbind(x,x[1,])
			}
			Neu <- rbind(Neu,x)
		}
		Slist[[i]] <- Neu
	}
	return(structure(Sou,SourceList=Slist,SAreas=SAreas,class=c("Sources","data.frame")))
}
