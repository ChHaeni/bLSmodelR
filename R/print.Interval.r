

print.Interval <- function(x,...){
	x[,12] <- sapply(x[,12],function(y){
		if(!is.na(y)){
			rplc <- sub(".{0,20}$","",sub("^.{20}","",y))
		} else {
			rplc <- "" 
		}
		if(rplc!=""){
			sub(rplc,"//...//",y,fixed=TRUE)
		} else {
			y
		}
	}, USE.NAMES = FALSE)
	x[,13] <- sapply(x[,13],function(y){
		if(!is.na(y)){
			rplc <- sub(".{0,20}$","",sub("^.{20}","",y))
		} else {
			rplc <- "" 
		}
		if(rplc!=""){
			sub(rplc,"//...//",y,fixed=TRUE)
		} else {
			y
		}
	}, USE.NAMES = FALSE)
	print.data.frame(x)
}


