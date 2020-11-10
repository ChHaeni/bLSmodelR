getCatalogs <- function(x,i=seq_len(nrow(x)),rn=NULL,Sensor=NULL){
  # convert old versions 
  sx <- as.character(substitute(x))
	x <- copy(x)
	setDT(x)
 	if(is.null(attr(x, "Version"))){
		warning(paste0("Object '", sx[min(length(sx), 2)], "' has not yet been converted to version 4.2+"))
		convert(x)
  }
 	switchNames(x)
	if(is.null(rn)&is.null(Sensor)){
		y <- x[i, .(rn, Sensor)]
	} else {
		rnI <- rn
		SensorI <- Sensor
		setkey(x,rn,Sensor)
		y <- x[.(rnI,SensorI), .(rn, Sensor)][1,]
	}
	setkey(y,rn,Sensor)
	out <- attr(x, "CalcSteps")[y,{
		CS <- unlist(strsplit(Calc.Sensor,","))
		.(Point.Sensor = CS,
			Catalog = paste(attr(x, "CatPath"), rep(Cat.Name, length(CS)), sep = "/"),
			seed = rep(seed, length(CS))
			)
	},by=.(rn, Sensor, Calc.Sensor)][,Calc.Sensor:=NULL]
	setnames(out,c("rn","Sensor","Calc.Sensor","Catalog","seed"))
	setkey(out,rn,Sensor,Calc.Sensor)		
	
	out
}
