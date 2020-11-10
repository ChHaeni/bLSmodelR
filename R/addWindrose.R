addWindrose <- function(WD,pos=2,frac=20,scF=1,nord=0,str.rot=FALSE){
	WD <- 90 - WD
	stR <- 28
	usr <- par()$usr
	dx <- diff(usr[1:2])
	dy <- diff(usr[3:4])
	R <- min(c(dx,dy))/frac
	switch(as.character(pos),
		"1"={corner <- usr[c(1,4)];Rx <- R;Ry <- -R;ddx <- 0.04*20/frac*dx;ddy <- -0.04*20/frac*dy},
		"3"={corner <- usr[c(2,3)];Rx <- -R;Ry <- R;ddx <- -0.04*20/frac*dx;ddy <- 0.04*20/frac*dy},
		"4"={corner <- usr[c(1,3)];Rx <- R;Ry <- R;ddx <- 0.04*20/frac*dx;ddy <- 0.04*20/frac*dy},
		{corner <- usr[c(2,4)];Rx <- -R;Ry <- -R;ddx <- -0.04*20/frac*dx;ddy <- -0.04*20/frac*dy}
		)
	M <- corner + c(Rx,Ry) + c(ddx,ddy)
	circ <- cbind(M[1] + cos(seq(0,2*pi,length.out=101)[-1])*R,M[2] + sin(seq(0,2*pi,length.out=101)[-1])*R)
	xEr <- xE <- M[1]+R*stR/frac;yEr <- yE <- M[2]
	xWr <- xW <- M[1]-R*stR/frac;yWr <- yW <- M[2]
	xNr <- xN <- M[1];yNr <- yN <- M[2]+R*stR/frac
	xSr <- xS <- M[1];ySr <- yS <- M[2]-R*stR/frac

	if(nord%%360!=0){
		nordRad <- nord/180*pi
		xEr <- M[1] + (xE-M[1])*cos(-nordRad) - (yE-M[2])*sin(-nordRad);yEr <- M[2] + (xE-M[1])*sin(-nordRad) + (yE-M[2])*cos(-nordRad)
		xWr <- M[1] + (xW-M[1])*cos(-nordRad) - (yW-M[2])*sin(-nordRad);yWr <- M[2] + (xW-M[1])*sin(-nordRad) + (yW-M[2])*cos(-nordRad)
		xNr <- M[1] + (xN-M[1])*cos(-nordRad) - (yN-M[2])*sin(-nordRad);yNr <- M[2] + (xN-M[1])*sin(-nordRad) + (yN-M[2])*cos(-nordRad)
		xSr <- M[1] + (xS-M[1])*cos(-nordRad) - (yS-M[2])*sin(-nordRad);ySr <- M[2] + (xS-M[1])*sin(-nordRad) + (yS-M[2])*cos(-nordRad)
	}

	# draw rose:
	lines(circ)
	srt=NA
	if(str.rot)srt <- -nord
	text(xEr,yEr,"E",cex=0.6*20/frac)
	text(xWr,yWr,"W",cex=0.6*20/frac)
	text(xNr,yNr,"N",cex=0.6*20/frac)
	text(xSr,ySr,"S",cex=0.6*20/frac)

	# arrow:
	ex1 <- 1.05*scF
	ex2 <- 1.15*scF
	dWD <- 5*scF
	ue <- c(x=M[1] + cos(WD/180*pi) * R*ex1,y=M[2] + sin(WD/180*pi) * R*ex1)
	ue1 <- c(x=M[1] + cos((WD+dWD)/180*pi) * R*ex2,y=M[2] + sin((WD+dWD)/180*pi) * R*ex2)
	ue2 <- c(x=M[1] + cos((WD-dWD)/180*pi) * R*ex2,y=M[2] + sin((WD-dWD)/180*pi) * R*ex2)
	le <- c(x=M[1] + cos(WD/180*pi-pi) * R*ex1,y=M[2] + sin(WD/180*pi-pi) * R*ex1)
	polygon(c(le[1],ue[1],ue1[1]),c(le[2],ue[2],ue1[2]),border=1)
	polygon(c(le[1],ue[1],ue2[1]),c(le[2],ue[2],ue2[2]),border=1,col=1)
}


