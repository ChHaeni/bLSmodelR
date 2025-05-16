addScaleBar <- function(pos=3,scale=1,units="meters",frac=5,cuts=c(1,2,5),minors=FALSE,cex=0.6){
	# scale: scale bar measures will be scaled by 1:scale (e.g using scale = 2 results in a measure of 1:2, i.e. 1.0cm -> display: 0.5cm)
	# frac: fraction of plot region x-length used for the length of scale
	usr <- par()$usr
	dx <- diff(usr[1:2])
	dy <- diff(usr[3:4])
	if(inherits(scale,"staticMap")){
		if(!requireNamespace("geosphere", quietly = TRUE)){
			stop("please install package geosphere: install.packages('geosphere')")
		}
		if(!requireNamespace("RgoogleMaps", quietly = TRUE)){
			stop("please install package RgoogleMaps: install.packages('RgoogleMaps')")
		}
		USR <- XY2LatLon(scale,usr[1:2],usr[3:4])
		scale <- sqrt((usr[1]-usr[2])^2+(usr[3]-usr[4])^2)/distGeo(USR[1,2:1],USR[2,2:1])
	}
	R <- numeric(3)
	R[1] <- round125(dx/frac/scale,cuts)
	R[2] <- round125(R[1]/2,cuts)
	R[3] <- round125(R[1]/5,cuts)
	R <- R*scale
	mR2 <- floor(R[1]/R[2]) - 1
	mR3 <- floor(R[1]/R[3]) - 1
	B <- R[1]/20/dx*dy
	switch(as.character(pos),
		"1"={corner <- usr[c(1,4)];R1 <- 0;R2 <- R[1];ddx <- 0.04*1.5*dx;ddy <- -0.04*1.5*dy},
		"2"={corner <- usr[c(2,4)];R1 <- -R[1];R2 <- 0;ddx <- -0.04*1.5*dx;ddy <- -0.04*1.5*dy},
		"3"={corner <- usr[c(2,3)];R1 <- -R[1];R2 <- 0;ddx <- -0.04*1.5*dx;ddy <- 0.04*1.5*dy},
		{corner <- usr[c(1,3)];R1 <- 0;R2 <- R[1];ddx <- 0.04*1.5*dx;ddy <- 0.04*1.5*dy}
		)
	lx <- corner[1] + ddx + R1
	rx <- corner[1] + ddx + R2
	ly <- corner[2] + ddy
	my <- ly + B
	polygon(c(lx,rx,rx,lx),c(ly,ly,my,my),col="white")
	lines(rep(lx,2),my+c(0,B/3*2))
	lines(rep(rx,2),my+c(0,B/3*2))
	rxH <- lx + R[2]
	polygon(c(rxH,rx,rx,rxH),c(ly,ly,my,my),col=1,border=1)
	rxS <- lx + R[3]
	polygon(c(lx,rxS,rxS,lx),c(ly,ly,my,my),col=1,border=1)
	if(!minors){
		mR3 <- mR2 <- 1
	}
	for(i in 1:mR2){
		lines(rep(lx + R[2]*i,2),my+c(0,B/3*2))
	}
	for(i in 1:mR3){
		lines(rep(lx + R[3]*i,2),my+c(0,B/3))
	}
	R <- R/scale
	text(lx,my+B+strheight("0")*cex,"0",cex=cex,adj=c(0.5,1))
	text(rxS,my+B+strheight("0")*cex,R[3],cex=cex,adj=c(0.5,1))
	text(rxH,my+B+strheight("0")*cex,R[2],cex=cex,adj=c(0.5,1))
	text(rx,my+B+strheight("0")*cex,R[1],cex=cex,adj=c(0.5,1))
	text((rx+lx)/2,my-B-strheight("0")*cex,units,cex=cex,adj=c(0.5,0))
}
