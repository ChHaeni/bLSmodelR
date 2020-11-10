tagInside <- function(Ctlg,Src,Sens=c("x-Coord (m)"=0,"y-Coord (m)"=0),tagBySourceName=tagPoly,tagPoly=FALSE){
	Ctlg[,tagInside:=if(tagBySourceName){rep("",nrow(Ctlg))}else{rep(FALSE,nrow(Ctlg))}]
	SourceAreaRelative <- data.table(Src)
	setnames(SourceAreaRelative,c("area","x","y","pid"))
	SourceAreaRelative[,":="(x=x-Sens[,"x-Coord (m)"],y=y-Sens[,"y-Coord (m)"])]
	tagNear(Ctlg,SourceAreaRelative)
	Ctlg[,inside0:=inside]
	Ctlg[,rn:=.I]
	setkey(Ctlg,rn)

	if(Ctlg[,any(inside)]){
		# tag Inside Source
		if(tagBySourceName){
			if(tagPoly){
				TDinside <- SourceAreaRelative[,
				{
					tagNear(Ctlg[,inside:=inside0],.(x=x,y=y))
					cbind(ID=Ctlg[(inside),rn],pnt.in.poly(Ctlg[(inside),cbind(x,y)],cbind(x,y)))
				},by=.(area,pid)][pip==1L,paste0(paste(area,pid,sep=","),collapse=";"),by=ID]	
			} else {
				TDinside <- SourceAreaRelative[,
				{
					tagNear(Ctlg[,inside:=inside0],.(x=x,y=y))
					cbind(ID=Ctlg[(inside),rn],pnt.in.poly(Ctlg[(inside),cbind(x,y)],cbind(x,y)))
				},by=.(area,pid)][pip==1L,paste0(area,collapse=","),by=ID]					
			}
			Ctlg[TDinside,tagInside:=V1]
		} else {
			TDinside <- SourceAreaRelative[,
			{
				tagNear(Ctlg[,inside:=inside0],.(x=x,y=y))
				cbind(ID=Ctlg[(inside),rn],pnt.in.poly(Ctlg[(inside),cbind(x,y)],cbind(x,y)))
			},by=.(area,pid)][,sum(pip),by=ID]
			Ctlg[TDinside,tagInside:=as.logical(V1)]
		}
	}
	Ctlg[,":="(inside=NULL,inside0=NULL,rn=NULL)]

	return(invisible(Ctlg))
}


