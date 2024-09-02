writeCatalog <- function(Ctlg, Name, compactTDcat = TRUE, ...) {
	# reduce sizes
	if (compactTDcat) compactCatalog(Ctlg, ...)

    # serialize first
    ## header (raw or serialized?)
    h_r <- attr(Ctlg, 'header')[[1]]
    h_s <- qs::qserialize(h_r, 'balanced')
    h_s_length <- length(h_s)
    ## entire catalog
    c_s <- qs::qserialize(Ctlg, 'balanced')
    c_s_length <- length(c_s)

    # write binary
    # open file connection
    on.exit(close(con))
    con <- file(Name, open = 'wb')
    # flag for new format
    writeBin('A', con, endian = 'little')
    # length of header
    writeBin(as.integer(h_s_length), size = 4, con, endian = 'little')
    # length of catalog
    writeBin(as.numeric(c_s_length), size = 8, con, endian = 'little')
    # header
    writeBin(h_s, con, endian = 'little')
    # catalog
    writeBin(c_s, con, endian = 'little')
    # TODO: update catfile
    # -> lock file (.00LOCK-catfile) -> same in rebuildCatListFile()
}
