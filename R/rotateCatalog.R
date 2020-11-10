rotateCatalog <- function(Cat,WD,back=FALSE){
	if(back){
		WDrad <- (180-WD)/180*pi
	} else {
		WDrad <- WD/180*pi
	}
	Cat[,":="(x=-x * sin(WDrad) + y*cos(WDrad),y=-x * cos(WDrad) - y*sin(WDrad))]
	invisible(Cat)
}
