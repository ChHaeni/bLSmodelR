extractResult <- function(DF,sortArgs,keep=NULL,quiet=FALSE,dropAttr=TRUE,na.last = FALSE){
	DFName <- deparse(substitute(DF))
	if(!inherits(DF,"bLSresult"))stop(DFName," does not inherit from class \"bLSresult\"!")
	if(!quiet)cat("~~~~~\nExtracting Results:\n")
	mc <- as.list(match.call()[["sortArgs"]][-1])
	atts <- attributes(DF)
	Names <- names(DF)
	if(!(isDT <- inherits(DF,"data.table"))){
		rn_remove_this_later_first <- row.names(DF)
		DF <- as.data.table(DF,keep.rownames = FALSE)
		DF[,rn_remove_this_later_first:=rn_remove_this_later_first]
		Names <- c(Names,"rn_remove_this_later_first")
	} else {
		DF <- copy(DF)
	}
		
	if(!quiet&length(mc)){
		sAn <- names(mc)
		if(is.null(sAn))sAn <- rep("",length(mc))
		empty <- sAn==""
		sAn[empty] <- as.character(mc[empty])
		dil <- rep("in increasing order",length(mc))
		dil[grepl("^[-]",sAn)] <- "in decreasing order"
		sAn <- gsub("^[-]","",sAn)
		sAn <- gsub("^[\"]","",sAn)
		sAn <- gsub("[\"]$","",sAn)
		dil[!empty] <- "by supplied list"
		cat("\nHierachically sorting",DFName,"by columns:\n")
		for(i in seq_along(sAn)){
			cat(paste0("\tlevel ",i,": Column ",paste0("\"",sAn[i],"\"")," of type \"",DF[,class(get(sAn[i]))],"\" sorted ",dil[i],"\n"))
			if(!empty[i]){
				kp <- eval(mc[[i]])
				rem <- DF[,unique(get(sAn[i])) %w/o% kp]
				if(length(rem))cat("\t\t- Removed column entry:",paste0(paste0("\"",rem,"\""),collapse=" ,"),"\n")
			}
		}		
	}

	if(length(mc))DF <- do.call("sortData",c(DF=quote(DF),mc,na.last=na.last))


	if(!is.null(keep)){
		nmspre <- c("rn_remove_this_later_first","Original Sonic Row","rn","Sensor","Source") %w/o% keep
		nms <- intersect(c(nmspre,keep),Names)
		DF <- DF[,nms,with=FALSE]
		if(!quiet){
			kOut <- paste0(DFName,"[,",paste0(paste0("\"",intersect(keep,Names),"\""),collapse=" ,"),"]")
			cat("\t-> Keeping Columns:",kOut,"\n\tdropped residual columns...\n")
		}
	}

	if(dropAttr){
		atts$class <- atts$class %w/o% "bLSresult"
		setattr(DF,"CalcSteps",NULL)
		setattr(DF,"CatPath",NULL)
		setattr(DF,"Catalogs",NULL)
		setattr(DF,"ModelInput", NULL)
		setattr(DF,"sessionInfo",NULL)
		setattr(DF,"ModelRunTime",NULL)
	} else {
		setattr(DF,"CalcSteps",atts$CalcSteps)
		setattr(DF,"CatPath",atts$CatPath)
		setattr(DF,"Catalogs",atts$Catalogs)
		setattr(DF,"ModelInput", atts$ModelInput)
		setattr(DF,"sessionInfo",atts$sessionInfo)
		setattr(DF,"ModelRunTime",atts$ModelRunTime)
	}

	if(!isDT){
		rn <- DF[,rn_remove_this_later_first]
		DF[,rn_remove_this_later_first:=NULL]
		setattr(DF, "row.names", rn)
		setattr(DF,"class",atts$class)
        setattr(DF, "sorted", NULL)
        setattr(DF, ".internal.selfref", NULL)		
	} else {
		setattr(DF,"class",atts$class)
    }

    if(!quiet)cat("~~~~~\n")

	return(DF)

}

