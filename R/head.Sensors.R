head.Sensors <- function(x,...){
  # convert old versions 
 	if(is.null(attr(x, "Version"))){
		warning("Object has not yet been converted to version 4.2+", call. = FALSE)
		x <- convert(x)
  }
	X <- x
	if(length(x) == 8){
		N_sensors <- table(X[, "Sensor Name"])
		Single_sensor <- names(N_sensors)[N_sensors == 1]
		ps <- X[, "Sensor Name"] %in% Single_sensor
		X[X[, "Distance between Point-Sensors (m)"] == 0, "Distance between Point-Sensors (m)"] <- "*"
		X[ps, "Distance between Point-Sensors (m)"] <- "-"
		X[, "Distance between Point-Sensors (m)"] <- factor(X[, "Distance between Point-Sensors (m)"])
		X[X[, "Number of Point-Sensors"] == 0, "Number of Point-Sensors"] <- "*"
		X[ps, "Number of Point-Sensors"] <- "-"
		X[, "Number of Point-Sensors"] <- factor(X[, "Number of Point-Sensors"])
		nPs <- sum(N_sensors == 1)
		nLs <- sum(N_sensors > 1)
		piheights <- sort(unique(x[ps, "Sensor Height (m)"]))
		pih <- length(piheights)
		if(pih > 0){
			brks <- ceiling(pih / 10)
			pNms <- character(brks)
			for(j in seq_len(brks)){
				if(j == brks && pih %% 10){
					pNms[j] <- paste(piheights[seq_len(pih %% 10) + (j - 1) * 10], collapse = ", ")
				} else {
					pNms[j] <- paste(piheights[seq_len(10) + (j - 1) * 10], collapse = ", ")
				}
			}			
		}
		LSnames <- sort(unique(x[!ps, "Sensor Name"]))
		cat("******\nSensors data.frame:\n******\n")
		cat(" Number of Sensors:", nPs + nLs, "\n")
		cat(" - Point Sensors:", nPs, "\n")
		if(pih > 0)cat(paste0(" \t", pih, " unique height", if(pih > 1)"s", ": ",
			paste(pNms, collapse = "\n \t                    "), " (m above ground)\n"))
		cat(" - Line Sensors:", nLs, "\n")
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
		cat("******\n\n")
	}
	class(X) <- "data.frame"
	head(X,...)
}
