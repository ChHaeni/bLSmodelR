		
cleanTemporary <- function(keep=NULL){
	files <- ls(.GlobalEnv,pattern="RemoveTemporary")
	files <- files %w/o% keep
	if(length(files)){
		cat("Removing temporary Catalogs...\n")
		eval(parse(text=paste0("file.remove(",paste(files,collapse=","),")")))
		rebuildCatListFile(dirname(files[1]))
		rm(list=files,envir=.GlobalEnv)
	}
}


