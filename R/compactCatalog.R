compactCatalog <- function(Ctlg,as.int=TRUE,round.values=as.int,tdvalues.round=c(1,2,2,4),uvw.round=3){
    if(attr(Ctlg,"is.compact")){
        if(attr(Ctlg,"is.int")&(!as.int|round.values)){
            tdvalues <- 10^attr(Ctlg,"tdvalues.round")
            uvw <- 10^attr(Ctlg,"uvw.round")
            Ctlg[,":="(
                    Time = as.numeric(Time/tdvalues[1])
                    ,x = as.numeric(x/tdvalues[2])
                    ,y = as.numeric(y/tdvalues[3])
                    ,wTD = as.numeric(wTD/tdvalues[4])
                )]
            setattr(Ctlg,"uvw0",sweep(attr(Ctlg,"uvw0"),2,uvw,"/"))
            setattr(Ctlg,"is.int",FALSE)
        }
    }
    if(round.values){
        tdvalues.round <- rep.int(tdvalues.round,4)[1:4]
        names(tdvalues.round) <- c("Time","x","y","wTD")
        uvw.round <- rep.int(uvw.round,3)[1:3]
        names(uvw.round) <- c("u0","v0","w0")
        if(attr(Ctlg,"is.compact")){
            tdvalues.round <- pmin(tdvalues.round,attr(Ctlg,"tdvalues.round"))
            uvw.round <- pmin(uvw.round,attr(Ctlg,"uvw.round"))
        }
        if(as.int){
            tdvalues.round <- pmin(tdvalues.round,c(4,6,6,8))
            tdvalues <- 10^tdvalues.round
            uvw.round <- pmin(uvw.round,7)
            uvw <- 10^uvw.round
            # Compact Catalog values:
            Ctlg[,":="(
                    Time = as.integer(round(Time,tdvalues.round[1])*tdvalues[1])
                    ,x = as.integer(round(x,tdvalues.round[2])*tdvalues[2])
                    ,y = as.integer(round(y,tdvalues.round[3])*tdvalues[3])
                    ,wTD = as.integer(round(wTD,tdvalues.round[4])*tdvalues[4])
                )]
            # uvw values:
            setattr(Ctlg,"uvw0",apply(sweep(round(attr(Ctlg,"uvw0"),rep(uvw.round,each=NROW(attr(Ctlg,"uvw0")))),2,uvw,"*"),2,as.integer))
            setattr(Ctlg,"is.int",TRUE)
        } else {
            # Compact Catalog values:
            Ctlg[,":="(
                    Time = round(Time,tdvalues.round[1])
                    ,x = round(x,tdvalues.round[2])
                    ,y = round(y,tdvalues.round[3])
                    ,wTD = round(wTD,tdvalues.round[4])
                )]
            # uvw values:
            setattr(Ctlg,"uvw0",round(attr(Ctlg,"uvw0"),rep(uvw.round,each=NROW(attr(Ctlg,"uvw0")))))
            setattr(Ctlg,"is.int",FALSE)

        }
        setattr(Ctlg,"tdvalues.round",tdvalues.round)
        setattr(Ctlg,"uvw.round",uvw.round)
        setattr(Ctlg,"is.compact",TRUE)
    } else if(as.int){
        stop("round values when coercing to integer...")
    }
    invisible(Ctlg)
}
check.compact <- function(Ctlg){
    out <- list(is.compact=FALSE,is.int=NULL,tdvalues.round=NULL,uvw.round=NULL)
    if(out$is.compact <- attr(Ctlg,"is.compact")){
        out$is.int <- attr(Ctlg,"is.int")
        out$tdvalues.round <- attr(Ctlg,"tdvalues.round")
        out$uvw.round <- attr(Ctlg,"uvw.round")
    }
    return(out)
}
