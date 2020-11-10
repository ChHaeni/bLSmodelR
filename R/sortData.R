
sortData <- function(DF,...,na.last =  FALSE){
	# (default) increasing, - descending, =c() 
	sortArgs_dont_use_that_in_colnames <- as.list(match.call(expand.dots=FALSE))[-1L][["..."]]
	sA_names <- names(sortArgs_dont_use_that_in_colnames)

	originalNames <- names(DF)
	if(isDT <- inherits(DF,"data.table")){
		DF <- copy(DF)
	} else {
		cl <- class(DF)
		rn_remove_this_later <- row.names(DF)
		DF <- as.data.table(DF,keep.rownames = FALSE)
		DF[,rn_remove_this_later:=rn_remove_this_later]
	}
	setattr(DF, "names", make.names(names(DF), unique = TRUE))
	newNames <- names(DF)

	if(!is.null(sA_names)){
		noname1 <- sA_names %in% ""
		noname <- which(noname1)
		nms_this_should_be_unique <- sA_names[!noname1]
		nms_this_should_be_unique_new <- newNames[chmatch(nms_this_should_be_unique,originalNames,nomatch=0)]
		for(index_that_is_unique in seq_along(nms_this_should_be_unique)){
			DF[,paste(nms_this_should_be_unique_new[index_that_is_unique],"ordering_and_remove_that_later",sep="_"):=chmatch(get(nms_this_should_be_unique_new[index_that_is_unique]),eval(sortArgs_dont_use_that_in_colnames[[nms_this_should_be_unique[index_that_is_unique]]]))]
			sortArgs_dont_use_that_in_colnames[[nms_this_should_be_unique[index_that_is_unique]]] <- paste(nms_this_should_be_unique_new[index_that_is_unique],"ordering_and_remove_that_later",sep="_")
		}
	} else {
		noname <- seq_len(length(sortArgs_dont_use_that_in_colnames))
	}

	sorting <- lapply(sortArgs_dont_use_that_in_colnames,function(x){
		if(is.call(x)){
			y2 <- -1L
			y <- as.character(x)[2]
		} else {
			y <- as.character(x)
			y2 <- 1L
		}
		list(y,y2)
	})

	scol <- sapply(sorting,"[[",1)
	
	for(i in noname)scol[i] <- newNames[chmatch(scol[i],originalNames,nomatch=NA)]
	if(any(is.na(scol)))stop("Some column names to sort data do not exist!")

	setorderv(DF,scol,sapply(sorting,"[[",2),na.last=na.last)

	if(!is.null(sA_names)){
		DF <- na.omit(DF,cols=paste(nms_this_should_be_unique_new,"ordering_and_remove_that_later",sep="_"))
		DF[,paste(nms_this_should_be_unique_new,"ordering_and_remove_that_later",sep="_"):=NULL]
	}

	if(!isDT){
		rn <- DF[,rn_remove_this_later]
		DF[,rn_remove_this_later:=NULL]
		setnames(DF,originalNames)
		setattr(DF, "row.names", rn)
		setattr(DF,"class",cl)
        setattr(DF, "sorted", NULL)
        setattr(DF, ".internal.selfref", NULL)
	} else {
		setnames(DF,originalNames)
	}

	return(DF)
}
