plot.Sources <- function(x,...)siteMap(x,...)
plot.Sensors <- function(x,...)siteMap(x,...)
plot.InputList <- function(x,...)siteMap(x,...)
plot.bLSresult <- function(x,...)siteMap(x,...)
siteMap <- function(x,y=NULL,xlab="x-Coord",ylab="y-Coord", panel.first = grid(),
	polygon.args=list(lwd=2,col=c("#107833","#23974A","#3EA962","#62C182","#93DBAC")),
	points.args=list(pch=20,cex=0.5,col=1), 
	lines.args=list(lty=3),
	sources.text.args=list(font=2,cex=0.5),
	sensors.text.args=list(font=2,cex=.5,pos=2),asp=1,add=FALSE,...){
	nPoly <- names(polygon.args)
	nPt <- names(points.args)
	nLS <- names(lines.args)
	polygon.args <- polygon.args[!(nPoly %in% c("x","y"))]
	points.args <- points.args[!(nPt %in% c("x","y"))]
	lines.args <- lines.args[!(nLS %in% c("x","y"))]
	ind <- !match(c("lwd","col"),nPoly,nomatch=0)
	polygon.args[c("lwd","col")[ind]] <- list(lwd=2,col=c("#107833","#23974A","#3EA962","#62C182","#93DBAC"))[ind]
	ind <- !match(c("pch","cex","col"),nPt,nomatch=0)
	points.args[c("pch","cex","col")[ind]] <- list(pch=20,cex=0.5,col=1)[ind]
	ind <- !match("lty",nLS,nomatch=0)
	lines.args["lty"[ind]] <- list(lty=3)[ind]
	# sources.text.args
	nSouT<- names(sources.text.args)
	ind <- !match(c("font","cex"),nSouT,nomatch=0)
	sources.text.args[c("font","cex")[ind]] <- list(font=2,cex=0.5)[ind]
	# sensors.text.args
	nSenT <- names(sensors.text.args)
	ind <- !match(c("font","cex","pos"),nSenT,nomatch=0)
	sensors.text.args[c("font","cex","pos")[ind]] <- list(font=2,cex=0.5,pos=2)[ind]
	SeNames <- sensors.text.args[["labels"]]
	nSenT <- names(sensors.text.args)
	ind <- match("labels",nSenT,nomatch=0)
	if(ind>0)sensors.text.args <- sensors.text.args[-ind]
	nSouT<- names(sources.text.args)
	nSenT <- names(sensors.text.args)
	
	if(inherits(x, "bLSresult")){
	  # convert old versions   
	 	if(is.null(attr(x, "Version"))){
			warning("Object has not yet been converted to version 4.2+", call. = FALSE)
			x <- convert(copy(x))
	  }
		x <- attr(x, "ModelInput")
	}
	if(inherits(x,"InputList")){
	  # convert old versions 
	 	if(is.null(attr(x$Sensors, "Version"))){
			warning("Sensors entry has not yet been converted to version 4.2+", call. = FALSE)
			Snsrs <- convert(x$Sensors)[, c("Sensor Name", "Sensor ID", "Node", "x-Coord (m)", "y-Coord (m)")]
	  } else {
	  	Snsrs <- x$Sensors[, c("Sensor Name", "Sensor ID", "Node", "x-Coord (m)", "y-Coord (m)")]
	  }
		Srcs <- procSources(x$Sources)
		Slist <- attr(Srcs,"SourceList")
	} else {
		Snsrs <- vector("list", 5)
		Slist <- numeric(0)
		Srcs <- matrix(numeric(0),ncol=3)	
	}
	xy <- list(x=x,y=y)
	cls1 <- match("Sensors",c(class(x)[1],class(y)[1]),nomatch=0)
	cls2 <- match("Sources",c(class(x)[1],class(y)[1]),nomatch=0)
	if(cls1>0){
	  # convert old versions 
	 	if(is.null(attr(xy[[cls1]], "Version"))){
			warning("Sensors object has not yet been converted to version 4.2+", call. = FALSE)
			Snsrs <- convert(xy[[cls1]])[, c("Sensor Name", "Sensor ID", "Node", "x-Coord (m)", "y-Coord (m)")]
	  } else {
	  	Snsrs <- xy[[cls1]][, c("Sensor Name", "Sensor ID", "Node", "x-Coord (m)", "y-Coord (m)")]
	  }		
	}
	if(cls2>0){
		Srcs <- procSources(xy[[cls2]])
		Slist <- attr(Srcs,"SourceList")
	}
	setDT(Snsrs)
	setnames(Snsrs, c("name", "id", "node", "x", "y"))
	rangeX <- range(c(Srcs[,2], Snsrs[, x]))
	rangeY <- range(c(Srcs[,3], Snsrs[, y]))
	if(!add){
		plot(rangeX,rangeY,type="n",panel.first=panel.first,asp=asp,xlab=xlab,ylab=ylab,...)
	}
	if(length(Slist)){

		colpos <- match("col",names(polygon.args),nomatch=0)
		PolyCol <- rep(polygon.args[["col"]],length(Slist))
		polygon.args <- polygon.args[-colpos]
		pCl <- length(PolyCol)
		reme <- numeric(0)
		xpos <- match("x",nSouT,nomatch=0)
		if(!xpos){
			tx <- sapply(Slist,function(x)mean(x[-1,1]))
		} else {
			tx <- rep(sources.text.args[["x"]],length(Slist))
			reme <- c(reme,xpos)
		}
		ypos <- match("y",nSouT,nomatch=0)
		if(!ypos){
			ty <- sapply(Slist,function(x)mean(x[-1,2]))
		} else {
			ty <- rep(sources.text.args[["y"]],length(Slist))
			reme <- c(reme,ypos)
		}
		nams <- match("labels",nSouT,nomatch=0)
		if(!nams){
			tlabels <- names(Slist)
		} else {
			tlabels <- sources.text.args[["labels"]]
			reme <- c(reme,nams)
		}
		cols <- match("col",nSouT,nomatch=0)
		if(!cols){
			tcol <- 1
		} else {
			tcol <- rep(sources.text.args[["col"]],length(Slist))
			reme <- c(reme,cols)
		}
		if(length(reme))sources.text.args <- sources.text.args[-reme]
	
		for(i in 1:length(Slist)){
			SA <- Slist[[i]]
			pindex <- SA[,3]
			PolyU <- unique(pindex)
			for(j in 1:length(PolyU)){
				# do.call(lines,c(list(x=SA[pindex==PolyU[j],1],y=SA[pindex==PolyU[j],2],col=PolyCol[(i-1)%%pCl+1]),polygon.args))
				do.call(polygon,c(list(x=SA[pindex==PolyU[j],1],y=SA[pindex==PolyU[j],2],col=PolyCol[(i-1)%%pCl+1]),polygon.args))
			}
			do.call(text,c(list(x=tx[i],y=ty[i],labels=tlabels[i],col=tcol[(i-1)%%pCl+1]),sources.text.args))
		}
	}

	if(nrow(Snsrs) > 0){
		Snsrs[, {
			if(.N > 1){
				do.call(lines, c(list(x = x, y = y), lines.args))
			}
		}, by = .(name, id)]


		a <- Snsrs[, factor(paste(x, y))]
		la <- levels(a)
		if(is.null(SeNames)){
			SeNames <- sapply(seq_along(la),function(i)paste(Snsrs[as.numeric(a) == i, name],collapse="\n"))
		} else {
			if(SeNames == ""){
				SeNames <- rep(SeNames,length(la))
			} else {
				SeNames0 <- SeNames
				SeNames <- sapply(seq_along(la),function(i)paste(SeNames0[as.numeric(a)==i],collapse="\n"))			
			}
		}
		colpos <- match("col",names(points.args),nomatch=0)
		Ptcol <- rep(points.args[["col"]],length(la))
		points.args <- points.args[-colpos]
		reme <- numeric(0)
		xpos <- match("x",nSenT,nomatch=0)
		if(!xpos){
			tx <- Snsrs[, x]
		} else {
			tx <- rep(sensors.text.args[["x"]], Snsrs[, .N])
			reme <- c(reme,xpos)
		}
		ypos <- match("y",nSenT,nomatch=0)
		if(!ypos){
			ty <- Snsrs[, y]
		} else {
			ty <- rep(sensors.text.args[["y"]], Snsrs[, .N])
			reme <- c(reme,ypos)
		}
		cols <- match("col",nSenT,nomatch=0)
		if(!cols){
			tcol <- Ptcol
		} else {
			tcol <- rep(sensors.text.args[["col"]],length(la))
			reme <- c(reme,cols)
		}
		if(length(reme))sensors.text.args <- sensors.text.args[-reme]
		for(i in seq_along(la)){
			ii <- which(as.numeric(a)==i)[1]
			do.call(points,c(list(x= Snsrs[ii, x], y = Snsrs[ii, y],col=Ptcol[i]),points.args))
			do.call(text,c(list(x=tx[ii],y=ty[ii],labels=SeNames[i],col=tcol[i]),sensors.text.args))
		}
	}
}
