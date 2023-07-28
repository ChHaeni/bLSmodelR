readCatalog <- function(Name, as.is = FALSE, header_only = FALSE) {
    # read from binary file
    on.exit(close(con))
    con <- file(Name, open = 'rb')

    # check flag
    sig <- readBin(con, 'character', endian = 'little')
    if (sig != 'A') {
        return(NULL)
    }
    # get length of header
    r_h_length <- readBin(con, 'integer', size = 4, endian = 'little')
    # get length of catalog
    r_c_length <- readBin(con, 'numeric', size = 8, endian = 'little')
    # read header
    r_h <- readBin(con, 'raw', n = r_h_length, endian = 'little')
    if (header_only) {
        # unserialize header
        header <- qs::qdeserialize(r_h, strict = TRUE)
        class(header) <- c('TDhead', 'character')
        return(header)
    }

    # read catalog
    r_c <- readBin(con, 'raw', n = r_c_length, endian = 'little')
    # unserialize catalog
    Ctlg <- alloc.col(qs::qdeserialize(r_c, strict = TRUE))

    # decompact sizes
	if (!as.is & attr(Ctlg, "is.int")) compactCatalog(Ctlg, as.int = FALSE)
    Ctlg
}
