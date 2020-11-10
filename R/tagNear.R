tagNear <- function(Cat,Poly){
	if(!("inside" %in% names(Cat)))Cat[,inside := TRUE]
	UseMethod("tagNear",Poly)
}
tagNear.matrix <- function(Cat,Poly){
	Cat[(inside),inside:=(
		x <= max(Poly[,1]) &
		x >= min(Poly[,1]) &
		y <= max(Poly[,2]) &
		y >= min(Poly[,2]))
	]
	return(Cat)
}
tagNear.data.table <- function(Cat,Poly){
	Cat[(inside),inside:=(
		x <= Poly[,max(x)] &
		x >= Poly[,min(x)] &
		y <= Poly[,max(y)] &
		y >= Poly[,min(y)])
	]
	return(Cat)
}
tagNear.list <- function(Cat,Poly){
	Cat[(inside),inside:=(
		x <= max(Poly[[1]]) &
		x >= min(Poly[[1]]) &
		y <= max(Poly[[2]]) &
		y >= min(Poly[[2]]))
	]
	return(Cat)
}
