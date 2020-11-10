rebuildCatListFile <- function(C.Path,File=character(0),fromScratch=FALSE){
	
	Existing <- dir(C.Path,pattern="Cat_Zm.*_[0-9]{14}$")
	ExistingFull <- paste(C.Path,Existing,sep="/")
	Catfile <- paste0(C.Path,"/.CatList")
	if(fromScratch){
		file.remove(Catfile)
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
			write.table(CatList,file=paste0(C.Path,"/.CatList0"),row.names=FALSE,col.names=TRUE)
			if(file.exists(Catfile))file.remove(Catfile)
			file.rename(paste0(C.Path,"/.CatList0"),Catfile)
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
			write.table(CatList,file=paste0(C.Path,"/.CatList0"),row.names=FALSE,col.names=TRUE)
			if(file.exists(Catfile))file.remove(Catfile)
			file.rename(paste0(C.Path,"/.CatList0"),Catfile)
		}
	} else {
		if(file.exists(Catfile))invisible(suppressWarnings(file.remove(Catfile)))
		CatList <- as.data.frame(c(list(a=character(0)),rep(list(a=numeric(0)),13)),stringsAsFactors=FALSE)
		colnames(CatList) <- c("Name","N0","ZSens","Ustar","L","Zo","Su_Ustar","Sv_Ustar","bw","C0","kv","A","alpha","MaxFetch")
	}
	invisible(CatList)
}

