rebuildCatListFile <- function(C.Path,File=character(0),fromScratch=FALSE){
	
	Existing <- dir(C.Path,pattern="Cat_Zm.*_[0-9]{14}$")
	ExistingFull <- paste(C.Path,Existing,sep="/")
	CatfileOrig <- paste0(C.Path,"/.CatList")
	Catfile <- tempfile(paste0('CatList', sample(1000, 1)))
	if(!fromScratch){
        file.copy(CatfileOrig, Catfile, overwrite = TRUE)
	} else {
		if (file.exists(Catfile)) file.remove(Catfile)
    }

	if(length(Existing)){
		if(!file.exists(Catfile)){
			CatList <- as.data.frame(c(list(a=character(0)),rep(list(a=numeric(0)),13)),stringsAsFactors=FALSE)
			colnames(CatList) <- c("Name","N0","ZSens","Ustar","L","Zo","Su_Ustar","Sv_Ustar","bw","C0","kv","A","alpha","MaxFetch")
		} else {
			CatList <- read.table(Catfile,header=TRUE,as.is=TRUE,colClasses=c("character",rep("numeric",13)))
			# check erroneous
			CatList <- CatList[grepl("Cat_Zm.*_[0-9]{14}$",CatList[,1]),]
			# remove duplicates
			CatList <- CatList[!(CatList[,1] %in% unique(CatList[duplicated(CatList[,1]),1])),]
		}		
		CatNames <- CatList[,1] %w/o% File
		if(!all(CatNames%in%Existing)){
			CatList <- CatList[CatNames%in%Existing,]
			write.table(CatList,file=Catfile,row.names=FALSE,col.names=TRUE)
		}
		if(!all(exCat <- Existing%in%CatNames)){
			CatAdd <- data.frame(matrix(NA,nrow=sum(!exCat),ncol=ncol(CatList)),stringsAsFactors=FALSE)
			colnames(CatAdd) <- colnames(CatList)
			for(i in seq_along(ExCat <- ExistingFull[!exCat])){
				Cat <- readCatalog(ExCat[i])
				Head <- unlist(strsplit(attr(Cat,"header"),"\n"))[-1]
				Whead <- matrix(as.numeric(gsub(".*[=] ","",Head)),nrow=1)
				CatAdd[i,-1] <- Whead
				CatAdd[i,1] <- basename(ExCat[i])
			}
			CatList <- rbind(CatList,CatAdd)
			write.table(CatList,file=Catfile,row.names=FALSE,col.names=TRUE)
		}
        file.remove(CatfileOrig)
        file.rename(Catfile, CatfileOrig)
        file.remove(Catfile)
	} else {
		if(file.exists(CatfileOrig))invisible(suppressWarnings(file.remove(CatfileOrig)))
		CatList <- as.data.frame(c(list(a=character(0)),rep(list(a=numeric(0)),13)),stringsAsFactors=FALSE)
		colnames(CatList) <- c("Name","N0","ZSens","Ustar","L","Zo","Su_Ustar","Sv_Ustar","bw","C0","kv","A","alpha","MaxFetch")
	}
	invisible(CatList)
}

