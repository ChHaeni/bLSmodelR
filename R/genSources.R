genSources <- function(...){
	ArgList <- list(...)
	if(length(ArgList)==0)stop("no input supplied\n")
	if(class(ArgList[[1]][[1]])=="list"){
		ArgList <- c(ArgList[[1]],ArgList[-1])
	}
	Out <- data.frame("SourceArea Name"=character(0),"x-Coord (m)"=numeric(0),"y-Coord (m)"=numeric(0),"Polygon ID"=numeric(0),stringsAsFactors=FALSE,check.names=FALSE)
	cnames <- names(Out)
	Lnams <- names(ArgList)
	if(is.null(Lnams))Lnams <- rep("",length(ArgList))
	Lind <- Lnams %in% ""
	Lnams[Lind] <- paste0("NoName",seq_along(which(Lind)))
	names(ArgList) <- Lnams
	Names <- unique(Lnams)
	isDF <- sapply(ArgList,is.data.frame)
	for(i in seq_along(Names)){
		nindex <- which(Names[i]==Lnams)
		if(length(nindex)==1&&isDF[i]&&length(ArgList[[nindex]])==4){
			app <- ArgList[[nindex]]
			if(any(is.na(app)))stop(paste0("Supplied argument (",Names[i],") contains NA values!\n"))
			app[,1] <- as.character(app[,1])
			names(app) <- cnames
			Out <- rbind(Out,app)
		} else {
			for(j in seq_along(nindex)){
				al <- ArgList[[nindex[j]]]
				if(any(is.na(al)))stop(paste0("Supplied argument (",Names[i],") contains NA values!\n"))
				if(class(al)=="data.frame"){
					appxy <- data.frame(al$x,al$y)
				} else if(length(al)==2){
					appxy <- data.frame(al$x,al$y)
				} else {
					appxy <- switch(tolower(substr(al[[1]],1,1)),
						"p" = data.frame(al$x,al$y),
						"c" = {
							M <- if(is.null(al$M)) c(0,0) else al$M
							R <- if(is.null(al$R)) 10 else al$R
							N <- if(is.null(al$N)) 100 else al$N
							cbind(
							M[1] + cos(seq(0,2*pi,length.out=N+1)[-1])*R
							,M[2] + sin(seq(0,2*pi,length.out=N+1)[-1])*R)},
						"r" = cbind(x=c(al$x1,al$x1,al$x2,al$x2),y=c(al$y1,al$y2,al$y2,al$y1))
					)				
				}
				app <- data.frame(Names[i],appxy,j,stringsAsFactors=FALSE)
				names(app) <- cnames
				Out <- rbind(Out,app)
			}
		}
	}

	return(structure(Out,class=c("Sources","data.frame")))
}

# subset Sources by name
'[.Sources' <- function(x, i, j, ...) {
    if (!missing(i) && is.character(i)) {
        i <- which(x[[1]] %in% i)
    }
    out <- `[.data.frame`(x, i, j, ...)
    if (ncol(out) != 4 && inherits(out, 'data.frame')) {
        class(out) <- 'data.frame'
    }
    out
}
'[<-.Sources' <- function(x, i, j, value) {
    if (!missing(i) && is.character(i)) {
        i <- which(x[[1]] %in% i)
    }
    out <- `[<-.data.frame`(x, i, j, value)
    if (ncol(out) != 4 && inherits(out, 'data.frame')) {
        class(out) <- 'data.frame'
    }
    out
}
