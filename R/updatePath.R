updatePath <- function(x,NewName){
	if(!dir.exists(NewName))warning("Path to touchdown catalogs (\"",NewName,"\") does not exist!")
	UseMethod("updatePath",x)
}
updatePath.data.frame <- function(x,NewName){
	name <- substitute(x)
	attr(x,"CatPath") <- NewName
	assign(as.character(name),x,parent.frame(),inherits=TRUE)
	return(invisible(x))
}
updatePath.data.table <- function(x,NewName){
	setattr(x,"CatPath",NewName)
	return(invisible(x))
}
