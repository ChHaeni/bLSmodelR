



plotFootprint <- function(x,SensorName,rn=NULL,MyMap=NULL,type=c("CE","wCE","uCE"),use.avg=FALSE,use.sym=FALSE,use.var=TRUE,wTDcutoff=NULL,origin=NULL,
	dx=2,dy=dx,breaks=function(x)quantile(c(0,max(x)),c(0.01,0.1,0.5,0.9)),xlim=c(-100,100),ylim=c(-100,100),add=FALSE,alpha=0.3,axs=c("r","i"),
	main=NULL,asp=1,fill=TRUE,sub=NULL,bg.col=NULL,addSource=TRUE,showMax=FALSE,showSensor=TRUE,dispSname=showSensor, N0 = NULL,
	lpos=NULL,showPerc=FALSE,leg.bg.col=grey(0.9),cpal=NULL,decPlaces=0,sigNums=1,addWR=FALSE,WRpos=2,WRfrac=20,WRscale=1,xy_transform=NULL,transformArgs=NULL,
	useSTRtree=TRUE,avoidGEOS = FALSE,addSB=FALSE,SBpos=3,StaticMapArgs=NULL,showLegend=TRUE){

	if(!requireNamespace("maptools")){
		stop("please install package maptools: install.packages('maptools')")
	}
	if(!requireNamespace("rgeos")){
		stop("please install package rgeos: install.packages('rgeos')")
	}
	# Notizen: 
	#	- Sources spaeter hinzufuegen
	#	- Methode fuer Catalog hinzufuegen
	#	- arguments besser organisieren
	#	- main auch bei MyMap
	#	- xlim/ylim xy_transform WGS84 etc... besser organisieren!!!
	#	- xy_transform --> coords_2_WGS84 oder so
	#	- TrueProj als Argument, oder besser: if add=TRUE, par()$usr -> MyMap$size???
	sx <- as.character(substitute(x))

	if(inherits(x, "footprint")){
		# breaks=rev(brks)
		xy_transform <- x$xy_transform
		transformArgs <- x$transformArgs
		SensorPosition <- x$SensorPosition
		WDmean <- x$WDmean
		xlimOriginal <- x$xlimOriginal
		ylimOriginal <- x$ylimOriginal
		xm <- x$xm
		ym <- x$ym
		xy <- x$xy
		if(!add){
			stop("only add=TRUE supported when providing a footprint object")
		}
	} else {
	  # convert old versions 
		x <- copy(x)
		setDT(x)
		switchNames(x)
	  if(is.null(attr(x, "Version"))){
			warning(paste0("Object '", sx[min(length(sx), 2)], "' has not yet been converted to version 4.2+"))
			convert(x)
	  }

		# checks:
		if(missing(SensorName)){
			stop("SensorName must be provided")
		} else {
			getThisSensor <- as.character(SensorName[1])
			if(!nrow(x[Sensor==getThisSensor]))stop("No results for Sensor: ",getThisSensor," availabe!")
		}
		stopifnot(type[1] %in% c("CE","wCE","uCE"))


		if(is.null(rn)){
			getThisRows <- x[Sensor==getThisSensor,unique(rn)]
		} else {
			getThisRows <- as.integer(rn)
		}
		getThisSource <- x[rn%in%getThisRows&Sensor==getThisSensor,Source[1]]
		index <- x[,which(rn%in%getThisRows&Sensor==getThisSensor&Source==getThisSource)]

		WD <- x[index,WD]
		if(is.null(N0)){
			N0 <- x[index[1],N0]
		} else {
			if(any(x[index,N0]<N0))stop("Choose N0 smaller or equal to original N0!")
		}
		
		Ustar <- x[index,Ustar] 
		Suu <- x[index,sUu]
		Svu <- x[index,sVu]
		bw <- x[index,bw]
		z <- x[index,SensorHeight]
		L <- x[index,L]
		Zo <- x[index,Zo]
		Source <- attr(procSources(attr(x,"ModelInput")$Sources),"SourceList")

		if(length(index)>1){
			WDmean <- (360 + atan2(sum(sin(WD/180*pi)*Ustar)/sum(Ustar),sum(cos(WD/180*pi)*Ustar)/sum(Ustar))/pi*180) %% 360
		} else {
			WDmean <- WD
		}

		Calc.Sensors <- procSensors(attr(x,"ModelInput")$Sensors)$"Calc.Sensors"
		SensorPosition <- Calc.Sensors[Calc.Sensors[, "Sensor Name"] %chin% SensorName, 
			c("Sensor Name", "x-Coord (m)", "y-Coord (m)", "Sensor ID", "Node")]
		SPR <- rotate(SensorPosition[,2:3],Angle=-WDmean)

		if(is.null(origin)){
			origin <- colMeans(SensorPosition[,2:3])
		}
		if(!add){
		if(is.null(main))main <- switch(type[1],
			"CE"="C-Footprint",
			"wCE"="w'C'-Footprint",
			"uCE"="u'C'-Footprint")
			if(is.null(MyMap)){
				plot(SensorPosition[,2],SensorPosition[,3],type="n",xlim=xlim + origin[1],ylim=ylim + origin[2],xlab="",ylab="",asp=asp,yaxt="n",xaxt="n",main=main,sub=sub)
				pr <- par()
				if(!is.null(bg.col))polygon(pr$usr[c(1,2,2,1)],pr$usr[c(3,3,4,4)],col=bg.col)
				axis(1,at=atx<-pretty(pr$usr[1:2]))
				axis(2,at=aty<-pretty(pr$usr[3:4]))
			} else {
				do.call(PlotOnStaticMap,c(list(MyMap=MyMap),StaticMapArgs))
				pr <- par()
				add <- TRUE
			}
		} else {
			pr <- par()
		}

		if(!is.null(MyMap)){
			if(!("TrueProj" %in% names(StaticMapArgs))||!StaticMapArgs$TrueProj){
				xy_transformOriginal <- xy_transform
				if(is.null(xy_transform)){
					xy_transform <- function(x,y,mymap=MyMap,...){
						LatLon2XY.centered(mymap,y,x)
					}
				} else {
					xy_transform <- function(x,y,mymap=MyMap,...){
						xy <- xy_transformOriginal(x,y,...)
						LatLon2XY.centered(mymap,xy$y,xy$x)
					}
				}
			}
		}
		if(addSource){
			for(i in Source){
				pg <- i
				if(!is.null(xy_transform)){
					pg[1:2] <- do.call(xy_transform,c(list(x=pg[,1],y=pg[,2]),transformArgs))
				}
				upg <- unique(pg[,3])
				for(j in upg){
					ind <- which(pg[,3]==j)
					polygon(pg[ind,1],pg[ind,2],col="#006837")
				}
			}
		}
		# transform x/y:
		if(!is.null(xy_transform)){
			SensorPosition[, 2:3] <- do.call(xy_transform,c(list(x=SensorPosition[,2],y=SensorPosition[,3]),transformArgs))
		}
		if(showSensor){
			sp_dt <- setnames(as.data.table(attr(x,"ModelInput")$Sensors[, 
				c("Sensor Name", "x-Coord (m)", "y-Coord (m)", "Sensor ID", "Node")]), 
				c("name", "x", "y", "id", "node"))[name == getThisSensor]

			# cat("hier stimmt's noch nicht falls google maps!\n")
			
			if(dispSname){
				if(sp_dt[, length(unique(id)) > 1]){
					sp_dt[, {
						ind <- order(node)
						lines(x[ind], y[ind], lty = 2, cex = 1.5)
						points(x[ind], y[ind], pch = 20, cex = 1.5)
						text(mean(x), mean(y), paste0(name[1], " (", .BY$id, ")"), cex = 0.6, pos = 4)
					}, by = id]
				} else {
					sp_dt[, {
						ind <- order(node)
						lines(x[ind], y[ind], lty = 2, cex = 1.5)
						points(x[ind], y[ind], pch = 20, cex = 1.5)
						text(mean(x), mean(y), name[1], cex = 0.6, pos = 4)
					}, by = id]
				}
			} else {
				sp_dt[, {
					ind <- order(node)
					lines(x[ind], y[ind], lty = 2, cex = 1.5)
					points(x[ind], y[ind], pch = 20, cex = 1.5)				
				}, by = id]				
			}
		}	
		xlim <- pr$usr[1:2]
		ylim <- pr$usr[3:4]

		xlimOriginal <- xlim
		ylimOriginal <- ylim
		
		if(!is.null(xy_transform)){
			parscale1 <- max(c(xlim[1],ylim[1]))
			parscale2 <- max(c(xlim[2],ylim[2]))
			lim1 <-  optim(c(-1,-1),function(x)sum((unlist(do.call(xy_transform,c(list(x=x[1],y=x[2]),transformArgs))) - c(xlim[1],ylim[1]))^2),control=list(reltol=sqrt(.Machine$double.eps)/parscale1))$par
			lim2 <-  optim(c(1,1),function(x)sum((unlist(do.call(xy_transform,c(list(x=x[1],y=x[2]),transformArgs))) - c(xlim[2],ylim[2]))^2),control=list(reltol=sqrt(.Machine$double.eps)/parscale2))$par
			xlim <- c(lim1[1],lim2[1])
			ylim <- c(lim1[2],lim2[2])
		}

		xylim <- cbind(x=c(xlim,rev(xlim)),y=rep(ylim,each=2))
		xylim <- rotate(xylim,Angle=-WDmean)
		xb <- seq(min(xylim[,1])-10*dx,max(xylim[,1])+10*dx,dx)
		yb <- seq(min(xylim[,2])-10*dy,max(xylim[,2])+10*dy,dy)
		xm <- xb[-1] - dx/2 
		ym <- yb[-1] - dy/2
		xy <- matrix(0,nrow=length(xm),ncol=length(ym))

		# loop over index
		for(j in seq_along(index)){
			cat(j,"/",length(index),"\n")
			Catalogs <- getCatalogs(x,index[j])
			CatNames <- Catalogs[,Catalog]
			Subset_seed <- Catalogs[, seed]
			dummy <- ""
			Ustardummy <- -9999
			# loop over point sensors
			for(i in seq_along(CatNames)){
				# i <- 1
				xytemp <- matrix(0,nrow=length(xm),ncol=length(ym))
				if(dummy!=CatNames[i]|x[index[j],Ustar!=Ustardummy]){
					dummy <- CatNames[i]
					Catalog <- readCatalog(CatNames[i])
					initializeCatalog(x[index[j],], Catalog=Catalog)
					uvw <- uvw0(Catalog)
					# subset?
					if(attr(Catalog,"N0")>N0){
						env <- globalenv()
					    oseed <- env$.Random.seed
					    set.seed(Subset_seed,kind="L'Ecuyer-CMRG")			
						takeSub <- sort(sample.int(attr(Catalog,"N0"),N0))
				        if (is.null(oseed)) {
				            rm(list = ".Random.seed", envir = env)
				        } else {
				            assign(".Random.seed", value = oseed, envir = env)
				        }
						attCat <- attributes(Catalog)
						indexNew <- 1:N0
						names(indexNew) <- takeSub
						Catalog <- Catalog[Traj_ID %in% takeSub,]
						Catalog[,":="(Traj_ID=indexNew[as.character(Traj_ID)])]
						for(ac in (names(attCat) %w/o% c("names","row.names",".internal.selfref","uvw0")))setattr(Catalog,ac,attCat[[ac]])
						uvw <- uvw[takeSub,]
						setattr(Catalog,"uvw0",uvw)
						HeaderSub <- paste(attr(Catalog,"header"),"\n*** Subset of ",N0," Trajectories ***\n\n",sep="")
						class(HeaderSub) <- c("TDhead","character")
						setattr(Catalog,"N0",N0)
						setattr(Catalog,"header",HeaderSub)
					}


					if(!is.null(wTDcutoff)) Catalog[wTD<wTDcutoff,wTD:=wTDcutoff]

					Catalog[,CE:=2/wTD/N0]

					if(use.avg)Catalog[,CE:=mean(CE)]

					if(type[1]!="CE"){					
						if(type[1]=="wCE"){
							w0 <- uvw0(Catalog)[,"w0"]
						} else {
							w0 <- uvw0(Catalog)[,"u0"] - calcU(Ustar, Zo, L, z)
						}
						Catalog[,CE:=w0[Traj_ID]*CE]	
					}

					if(use.sym){
						Catalog <- Catalog[rep(1:.N,2)]
						Catalog[1:(.N/2),y:=-y]
						Catalog[,CE:=CE/2]
					}
					
				}
				Ctlg <- copy(Catalog)
				rotateCatalog(Ctlg,WD[j]-WDmean)
				Ctlg <- Ctlg[,":="(x=x+SPR[i,1],y=y+SPR[i,2],CE=CE/length(CatNames)/length(index))][
					x>=xb[1]&x<=xb[length(xb)]&y>=yb[1]&y<=yb[length(yb)]
				]

				Ctlg[,":="(xtag=findInterval(x,xb, rightmost.closed = TRUE),
					ytag=findInterval(y,yb, rightmost.closed = TRUE))]
				Cat <- Ctlg[,sum(CE),by=.(xtag,ytag)]
				xytemp[Cat[,cbind(xtag,ytag)]] <- Cat[,V1]

				if(use.var){
					xb2l <- Ctlg[,c(rev(seq(SPR[i,1]-3*dx/2,min(x,SPR[i,1])-10*dx,-2*dx)),seq(SPR[i,1]+dx/2,max(x,SPR[i,1])+10*dx,2*dx))]
					yb2l <- Ctlg[,c(rev(seq(SPR[i,2]-3*dy/2,min(y,SPR[i,2])-10*dy,-2*dy)),seq(SPR[i,2]+dy/2,max(y,SPR[i,2])+10*dy,2*dy))]

					xb2r <- xb2l + dx
					yb2r <- yb2l + dx
					Ctlg[,":="(
						xtag_ll=findInterval(x,xb2l),ytag_ll=findInterval(y,yb2l)
						,xtag_lr=findInterval(x,xb2l),ytag_lr=findInterval(y,yb2r)
						,xtag_rl=findInterval(x,xb2r),ytag_rl=findInterval(y,yb2l)
						,xtag_rr=findInterval(x,xb2r),ytag_rr=findInterval(y,yb2r)
						)]
					Cat_ll <- Ctlg[,sum(CE),by=.(xtag_ll,ytag_ll)]
					Cat_lr <- Ctlg[,sum(CE),by=.(xtag_lr,ytag_lr)]
					Cat_rl <- Ctlg[,sum(CE),by=.(xtag_rl,ytag_rl)]
					Cat_rr <- Ctlg[,sum(CE),by=.(xtag_rr,ytag_rr)]

					xy_ll <- xy_lr <- xy_rl <- xy_rr <- matrix(0,nrow=length(xb2r)-1,ncol=length(yb2r)-1)
					xy_ll[Cat_ll[,cbind(xtag_ll,ytag_ll)]] <- Cat_ll[,V1]
					xy_lr[Cat_lr[,cbind(xtag_lr,ytag_lr)]] <- Cat_lr[,V1]
					xy_rl[Cat_rl[,cbind(xtag_rl,ytag_rl)]] <- Cat_rl[,V1]
					xy_rr[Cat_rr[,cbind(xtag_rr,ytag_rr)]] <- Cat_rr[,V1]

					xmil <- findInterval(xm,xb2l,all.inside = TRUE)
					ymil <- findInterval(ym,yb2l,all.inside = TRUE)
					xmir <- findInterval(xm,xb2r,all.inside = TRUE)
					ymir <- findInterval(ym,yb2r,all.inside = TRUE)

					xyi_ll <- as.matrix(expand.grid(x=xmil,y=ymil,KEEP.OUT.ATTRS=FALSE))
					xyi_rr <- as.matrix(expand.grid(x=xmir,y=ymir,KEEP.OUT.ATTRS=FALSE))
					xyi_lr <- as.matrix(expand.grid(x=xmil,y=ymir,KEEP.OUT.ATTRS=FALSE))
					xyi_rl <- as.matrix(expand.grid(x=xmir,y=ymil,KEEP.OUT.ATTRS=FALSE))

					Mu <- (xy_ll[xyi_ll] + xy_lr[xyi_lr] + xy_rl[xyi_rl] + xy_rr[xyi_rr])/4
					Var <- (4*xytemp - Mu)^2 
					xytemp <- (xytemp*Var + (max(Var)-Var)/4*Mu)/(max(Var))
				}

				xy <- xy + xytemp/(dx*dy)

			}
		}
		rm(Catalog,Ctlg,Cat)
		xy[,c(1,ncol(xy))] <- 0
		xy[c(1,nrow(xy)),] <- 0
	}

	brks <- breaks(xy)

	cl <- contourLines(xm,ym,xy,levels=brks)
	clRot <- lapply(cl,rotate,Angle=WDmean)
	clFac <- as.factor(sapply(clRot,"[[","level"))
	uclFac <- unique(clFac)

	# transform x/y:
	if(!is.null(xy_transform)){
		clRot <- lapply(clRot,function(x){
			a <- do.call(xy_transform,c(list(x=x$x,y=x$y),transformArgs))
			list(level=x$level,x=a[[1]],y=a[[2]])
		})
	}

	clP <- lapply(clRot,function(x)sp::Polygon(cbind(x$x,x$y)))
	clPs <- lapply(uclFac,function(x,y,z,a,b)maptools::checkPolygonsHoles(sp::Polygons(y[z==x],as.character(x)),useSTRtree=a,avoidGEOS=b),y=clP,z=clFac,a=useSTRtree,b=avoidGEOS)
	clSp <- sp::SpatialPolygons(clPs,as.integer(uclFac))
	
	# if(is.null(xlim)&is.null(ylim)){
	# 	xyr <- rotate(cbind(x=range(xm),y=range(ym)),Angle=WDmean+180)
	# 	xlim <- range(xyr[,1]) + range(SensorPosition[,2])
	# 	ylim <- range(xyr[,2]) + range(SensorPosition[,3])
	# }

	# Loesung:
	mlt <- switch(axs[1],
		"i"=1.08,
		1
		)
	xlim2 <- xlimOriginal - diff(xlimOriginal)*(1-1/mlt)/2;ylim2 <- ylimOriginal - diff(ylimOriginal)*(1-1/mlt)/2
	c1 <- sp::Polygon(cbind(c(xlim2[1],xlim2[2],xlim2[2],xlim2[1], xlim2[1]),c(ylim2[1],ylim2[1],ylim2[2],ylim2[2],ylim2[1])))
	c2 <- sp::Polygons(list(c1), "sclip")
	Pclip <- sp::SpatialPolygons(list(c2))

	SpPclip <- rgeos::gIntersection(clSp, Pclip, byid = TRUE)
	
	alphachar <- as.hexmode(round(alpha*255))
	
	if(is.null(cpal)){
		cpal <- if(type[1]=="CE") ConcPalette(length(brks)) else FluxPalette(length(brks))
	} else {
		cpal <- cpal(length(brks))
	}

	if(!add){abline(h=aty,lty="dotted",col="lightgrey");abline(v=atx,lty="dotted",col="lightgrey")}

	if(fill)fillcol <- paste0(cpal,alphachar) else fillcol <- NA

	sp::plot(SpPclip,col=fillcol,border=cpal,add=TRUE)
	box()

    mi <- which(xy==max(xy),arr.ind=T)
	if(showMax){
		maxli <- rotate(cbind(x=xm[mi[1]],y=ym[mi[2]]),Angle=WD+90,Center=origin)+c(SensorPosition[,2],SensorPosition[,3])
		points(maxli,pch=20,cex=0.4,col="red")
		text(maxli[1],maxli[2],"max.",cex=0.6,pos=4,offset=0.2,col="red")
	}

	if(showLegend){
		if(is.null(lpos))lpos <- c("bottomleft","topleft","topright","bottomright")[floor(WDmean/90)+1]
		if(showPerc){
			ltex <- sprintf(paste0("%2.",decPlaces,"f%%"),rev(brks)/max(xy)*100)
			ltex <- paste0(">",ltex)
			temp <- legend(lpos,legend=rep(" ",length(ltex)),title="% of max.",text.width=strwidth(paste0("00",if(decPlaces>0){paste0(".",paste0(rep("0",decPlaces),collapse=""))}else{""}," %")),fill=paste0(rev(cpal),alphachar),border=rev(cpal),bg=leg.bg.col)
			text(temp$rect$left + temp$rect$w, temp$text$y,ltex, pos = 2)
		} else {
			ltex <- signif(rev(brks),sigNums)
			ltex <- paste0(">",ltex)
			legend(lpos,legend=ltex,xjust=0,fill=paste0(rev(cpal),alphachar),border=rev(cpal),bg=leg.bg.col)
		}
	}
	if(addWR)addWindrose(WDmean,pos=WRpos,frac=WRfrac,scF=WRscale)
	scale <- if(is.null(MyMap)) 1 else MyMap
	if(addSB)addScaleBar(SBpos,scale)
	wm <- which(xy==max(xy),arr.ind=TRUE)
	WDrad <- WDmean/180*pi
	out <- structure(list(
		breaks=rev(brks),
		range.z=range(xy),
		max.index = mi,
		xy_transform = xy_transform,
		transformArgs = transformArgs,
		SensorPosition = SensorPosition,
		WDmean = WDmean,
		xlimOriginal = xlimOriginal,
		ylimOriginal = ylimOriginal,
		xm = xm,
		ym = ym,
		xy = xy), class = "footprint"
		)
	return(invisible(out))

}

