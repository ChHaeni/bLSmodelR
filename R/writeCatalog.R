writeCatalog <- function(Ctlg, Name, compactTDcat = TRUE, ...) {
	# write data:
	if (compactTDcat) compactCatalog(Ctlg, ...)
	qs::qsave(Ctlg, Name, 'balanced')
}
