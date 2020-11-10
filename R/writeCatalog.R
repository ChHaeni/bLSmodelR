writeCatalog <- function(Ctlg,Name,compactTDcat=TRUE,...){
	# write data:
	if(compactTDcat)compactCatalog(Ctlg,...)
	saveRDS(Ctlg,Name)
}
