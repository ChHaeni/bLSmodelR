assignCat <- function(tdcat,tdlist){
	name = substitute(tdcat)
	attCat <- attributes(tdcat)
	tdcat <- rbindlist(list(tdcat,tdlist))
	for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref")))setattr(tdcat,ac,attCat[[ac]])
	assign(as.character(name), tdcat, parent.frame(), inherits = TRUE)
	invisible(tdcat)
}

