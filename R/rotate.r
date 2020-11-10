

rotate <- function(x,xAxisTo=0,Center=c(0,0),Angle=xAxisTo-90,nms=c("x","y")){
	ang <- -Angle/180
	ctr <- Center
	UseMethod("rotate",x)
}
rotate.default <- function(x,xAxisTo=0,Center=c(0,0),Angle=xAxisTo-90,nms=c("x","y")){
	if(!is.numeric(nms)){
		nms <- match(nms,names(x),nomatch=0)
		nms[!nms] <- (c(1,2) %w/o% nms)[seq_len(sum(nms==0))]		
	}
	x[[nms[1]]] <- x[[nms[1]]] - ctr[1]
	x[[nms[2]]] <- x[[nms[2]]] - ctr[2]
	out <- x
	out[[nms[1]]] <- x[[nms[1]]] * cospi(ang) - x[[nms[2]]] * sinpi(ang)  + ctr[1]
	out[[nms[2]]] <- x[[nms[1]]] * sinpi(ang) + x[[nms[2]]] * cospi(ang)  + ctr[2]
	return(out)
}
rotate.matrix <- function(x,xAxisTo=0,Center=c(0,0),Angle=xAxisTo-90,nms=c("x","y")){
	if(!is.numeric(nms)){
		nms <- match(nms,colnames(x),nomatch=0)
		nms[!nms] <- (c(1,2) %w/o% nms)[seq_len(sum(nms==0))]		
	}
	x[,nms[1]] <- x[,nms[1]] - ctr[1]
	x[,nms[2]] <- x[,nms[2]] - ctr[2]
	out <- x
	out[,nms[1]] <- x[,nms[1]] * cospi(ang) - x[,nms[2]] * sinpi(ang)  + ctr[1]
	out[,nms[2]] <- x[,nms[1]] * sinpi(ang) + x[,nms[2]] * cospi(ang)  + ctr[2]
	return(out)
}
rotate.data.table <- function(x,xAxisTo=0,Center=c(0,0),Angle=xAxisTo-90,nms=c("x","y")){
	if(is.numeric(nms)){
		nms <- names(x)[nms]
	}
	xnms <- x[,.(
		nms1 = get(nms[1]) - ctr[1]
		,nms2 = get(nms[2]) - ctr[2]
		)]
	x[,nms[1] := xnms[,nms1] * cospi(ang) - xnms[,nms2] * sinpi(ang)  + ctr[1],with=FALSE]
	x[,nms[2] := xnms[,nms1] * sinpi(ang) + xnms[,nms2] * cospi(ang)  + ctr[2],with=FALSE]
	invisible(x)
}

