join <- function(...)UseMethod("join")
join.Sensors <- function(...){
    allargs <- list(...)
    allargs <- allargs[sapply(allargs, inherits, what = "Sensors")]
    Nams <- unlist(lNams <- lapply(allargs, function(x) unique(x[, "Sensor Name"])))
    if(any(dNams <- duplicated(Nams))){
        warning("You supplied sensors with identical names. Duplicated sensor names will be treated as one (separated) sensor path!\n")
        uNams <- unique(Nams[dNams])
        for(i in uNams){
            index <- which(sapply(lNams, function(x) i %in% x))
            # extend ids
            for(j in seq_along(index)){
                temp <- allargs[[index[j]]]
                allargs[[index[j]]][temp[, "Sensor Name"] %in% i, "Sensor ID"] <- paste(
                    temp[temp[, "Sensor Name"] %in% i, "Sensor ID"], j, 
                    sep = "_")
            }
        }
    }
    structure(
    	do.call(rbind.data.frame, allargs)
        , class = c("Sensors", "data.frame")
        )
}
join.Sources <- function(...){
    allargs <- list(...)
    # allargs <- allargs[lengths(allargs) > 0L]
    allargs <- allargs[sapply(allargs,inherits,what="Sources")]
    Nams <- unlist(lNams<-lapply(allargs,function(x)unique(x[,1])))
    if(any(dNams <- duplicated(Nams))){
        warning("You supplied sources with identical names. Duplicated source names will be treated as one source!\n")
        uNams <- unique(Nams[dNams])
        for(i in uNams){
            index <- which(sapply(lNams,function(x)i %in% x))
            maxis <- cumsum(c(0,sapply(allargs[index],function(x)max(x[,4]))[-length(index)]))
            allargs[index] <- mapply(function(x,y){x[,4] <- x[,4]+y;return(x)},x=allargs[index],y=maxis,SIMPLIFY=FALSE)
        }
    }
    Out <- structure(
        do.call(rbind.data.frame,allargs)
        ,class=c("Sources","data.frame")
        )
    return(Out)
}
join.bLSresult <- function(...,asDT=TRUE){
    allargs <- list(...)
    allargs <- allargs[sapply(allargs,inherits,what="bLSresult")]
    checkDF <- !sapply(allargs,is.data.table)
    if(any(checkDF)){
        for(i in which(checkDF)){
            x <- copy(allargs[[i]])
            setDT(x)
            allargs[[i]] <- x
        }
    }
    # change rn and convert to 4.2+
    for(i in seq_along(allargs)){
        allargs[[i]][, rn := paste0(i, "_", rn)]
        attr(allargs[[i]], "Catalogs")[, rn := paste0(i, "_", rn)]
        attr(allargs[[i]], "CalcSteps")[, rn := paste0(i, "_", rn)]
        row.names(attr(allargs[[i]], "ModelInput")$"Interval") <- paste0(i, "_", row.names(attr(allargs[[i]], "ModelInput")$"Interval"))
        allargs[[i]] <- convert(allargs[[i]])
    }
    out <- rbindlist(allargs,fill = TRUE)
    setattr(out, "Version", "4.2+")
    # attributes:c("CalcSteps","CatPath","Catalogs","ModelInput","ModelRunTime")
    setattr(out,"ModelRunTime",Reduce("+",lapply(allargs,attr,which="ModelRunTime")))
    ### CatPath
    CPaths <- unique(sapply(allargs,attr,which="CatPath"))
    if(length(CPaths) > 1){
        CPaths <- CPaths[1]
        warning("bLS results do not share the same path to the catalog folder! Move your catalogs to:\n\t",CPaths)
    }
    setattr(out,"CatPath",CPaths)
    ### ModelInput
    InLists <- lapply(allargs,attr,which="ModelInput")
    MInput <- attr(allargs[[1]],"ModelInput")
    Sens <- unique(rbindlist(lapply(InLists,"[[","Sensors")))
    setnames(Sens, c("name", "id", "node", "x", "y", "z", "d", "n"))
    MInput$Sensors <- genSensors(Sens)
    # bind Sources & check for duplicated sources
    Sour <- as.data.frame(unique(rbindlist(lapply(InLists, function(x) {
                cbind(x[["Sources"]], row = seq_len(nrow(x[["Sources"]])))
    }))))
    # gen new unique sources
    MInput$Sources <- genSources(Sour[, 1:4])
    if(!all(sapply(lapply(InLists,"[[","Model"),function(x)identical(x,MInput$Model)))){
        warning("ignored differing (global) model parameters!")
    }
    if(!all(sapply(lapply(InLists,"[[","Tolerances"),function(x)identical(x,MInput$Tolerances)))){
        warning("ignored differing tolerance parameters!")
    }
    MInput$Interval <- as.data.frame(rbindlist(lapply(InLists,"[[","Interval"),fill = TRUE))
    setattr(out,"ModelInput",MInput)
    ### CalcSteps
    setattr(out,"CalcSteps",rbindlist(lapply(allargs,attr,which="CalcSteps"),fill = TRUE))
    ### Catalogs
    setattr(out,"Catalogs",rbindlist(lapply(allargs,attr,which="Catalogs"),fill = TRUE))
    # set Catalogs and CalcSteps keys
    setkey(attr(out, 'CalcSteps'),rn,Sensor)
    setkey(attr(out, 'Catalogs'),rn,Sensor,PointSensor)
    ### class bLSresult
    setattr(out, "class", c("bLSresult", "data.table", "data.frame"))
    if(asDT){
        out
    } else {
        setDF(out)
    }   
}

join.deposition <- function(...){
    allargs <- list(...)
    allargs <- allargs[sapply(allargs,inherits,what="deposition")]
    out <- rbindlist(allargs,fill = TRUE)
    ### vDep
    vdList <- lapply(allargs,attr,which="vDep")
    vDep <- sapply(vdList,"[[","vDep")
    names(vDep) <- NULL
    vDepSpatial <- lapply(vdList,"[[","vDepSpatial")
    setattr(out,"vDep",list(vDep=vDep,vDepSpatial=vDepSpatial))
    ### class bLSresult, deposition
    setattr(out, "class", c("deposition", "bLSresult", "data.table", "data.frame"))
    out
}

