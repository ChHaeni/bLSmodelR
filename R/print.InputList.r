

print.InputList <- function(x,...){
	X <- x
	cat("~~~~~~ Model Input List ~~~~~~\n\n")
	# Sensors:
	x <- X$Sensors
 	if(is.null(attr(x, "Version"))){
		warning("Entry 'Sensors' has not yet been converted to version 4.2+", call. = FALSE)
		x <- convert(x)
  }
	if(length(x) == 8){
		N_sensors <- table(x[, "Sensor Name"])
		Single_sensor <- names(N_sensors)[N_sensors == 1]
		ps <- x[, "Sensor Name"] %in% Single_sensor
		nPs <- sum(N_sensors == 1)
		xPs <- x[ps,]
		nLs <- sum(N_sensors > 1)
		piheights <- sort(unique(x[ps, "Sensor Height (m)"]))
		pih <- length(piheights)
		Pnames <- sapply(piheights, function(y){
			out <- paste(xPs[xPs[, "Sensor Height (m)"] == y, "Sensor Name"],collapse=",")
			rplc <- sub(".{0,30}$","",sub("^.{30}","",out))
			if(rplc!=""){
				sub(rplc,"//...//",out,fixed=TRUE)
			} else {
				out
			}
		})
		LSnames <- sort(unique(x[!ps, "Sensor Name"]))
		cat("*** Number of Sensors:", nPs + nLs, "\n")
		cat(" Point Sensors:", nPs, "\n")
		if(pih>0)cat(paste0("   ",pih," unique height",if(pih>1)"s",": \n   - ",paste(paste0(Pnames,": ",piheights," (m above ground)"),collapse="\n   - "),"\n"))
		cat(" Line Sensors:",nLs,"\n")
		if(nLs > 0){
			for(i in seq(LSnames)){
                ind <- match(LSnames[i], x[, 'Sensor Name'])
				l <- by(x[ind,], x[ind, "Sensor ID"], function(y){
					ord <- order(y[, "Node"])
					sum(sqrt(diff(y[ord, "x-Coord (m)"]) ^ 2 + diff(y[ord, "y-Coord (m)"]) ^ 2 + diff(y[ord, "Sensor Height (m)"]) ^ 2))
				})
				ds <- x[ind, "Distance between Point-Sensors (m)"]
				if(any(ds > 0)){
					dr <- unique(round(range(ds[ds > 0]), 1))
				} else {
					dr <- 0
				}
				cat(paste0(" \t ", LSnames[i], ":  ", paste(sort(range(x[ind, "Sensor Height (m)"])), collapse = " to ")
					," (m above ground)\t",round(sum(l), 1)," m (",sum(x[ind, "Number of Point-Sensors"])," PS. path: ",sum(x[ind, "Number of Point-Sensors"]) - 
					length(unique(x[ind, "Sensor ID"]))," x ",
					paste(dr, collapse = " to ")," m)\n"))
			}
		}
		cat("\n")
	} else {
		cat("xxxxxx\nInvalid Sensors data.frame!!!\nxxxxxx\n")
	}
	# Sources:
	x <- X$Sources
	if(length(x)==4){
		x <- procSources(x)
		uS <- unique(x[,1])
		nS <- length(uS)
		cat("*** Number of Sources:",nS,"\n")
		SourceList <- attr(x,"SourceList")
		for(i in seq(nS)){
			p <- SourceList[[uS[i]]]
			np <- length(unique(p[,3]))
			cat("   -",uS[i],":",round(sum(by(p,p[,3],getArea)),2),"(m2)",if(np>1){paste0("(",length(unique(p[,3]))," polygons)")},"\n")
		}
		cat("\n")
	} else if(is.null(x)){
		cat("*** Number of Sources: 0\n")
		cat("No Sources defined!\n\n")
	} else {
		cat("xxxxxx\nInvalid Sources data.frame!!!\nxxxxxx\n")
	}
	# Model:
	x <- X$Model
	cat("*** Model Parameters:\n")
	for(i in names(x))cat("    ",i,"= ",x[[i]],"\n")
	cat("\n")
	# Tolerances:
	x <- X$Tolerances
	rn <- rownames(x)
	cat("*** Tolerances (in %):\n")
	for(i in 1:6)cat("    ",rn[i],"= ",x[i,1],"\n")
		cat("\n")
	# Interval:
	x <- X$Interval
	cat("*** Interval Input:\n")
    print(x, show_all = FALSE, n_split = 3)
}


