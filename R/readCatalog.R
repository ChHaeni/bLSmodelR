readCatalog <- function(Name,as.is=FALSE){
	# get data:
	Ctlg <- alloc.col(readRDS(Name))
	if(!as.is&attr(Ctlg,"is.int"))compactCatalog(Ctlg,as.int=FALSE)
	return(Ctlg)
}
