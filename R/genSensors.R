genSensors <- function(...){
	ArgList <- list(...)
	if (length(ArgList) == 0) stop("no input supplied\n")
	Out <- NULL
	Names <- names(ArgList)
	isDF <- sapply(ArgList, is.data.frame)
	if(is.null(Names)) Names <- rep("", length(ArgList))
	isNull <- Names %in% ""
	Names[isNull] <- as.character(which(isNull))
	names(ArgList) <- Names
	names(isNull) <- Names
	DFNames <- unlist(lapply(ArgList[isDF], function(x) {
		if(is.character(x[1, 1]) || is.factor(x[1, 1])){
			unique(as.character(x[,1]))
		} else if(any(Sn <- (c("name", "sensor", "Sensor Name", "sensor_name") %in% names(x)))){
			unique(as.character(x[, c("name", "sensor", "Sensor Name", "sensor_name")[Sn]]))
		} else {
			NULL
		}
	}))
	if(any(ResName <- grepl(".*[.]{1}[0-9]*$", c(Names[!isDF], DFNames)))){
			if(sum(ResName)>1){
				stop("Sensor Names are not allowed to end with .[0-9]*!\n\t -> Sensor Names: ",paste(unique(c(Names[!isDF],DFNames)[ResName]),collapse=", ")," are not valid!\n")
			} else {
				stop("Sensor Names are not allowed to end with .[0-9]*!\n\t -> Sensor Name: ",paste(unique(c(Names[!isDF],DFNames)[ResName]),collapse=", ")," is not valid!\n")
			}
	}
	if(Dups <- any(duplicated(c(Names[!isDF],DFNames)))){
		dupNames <- unique(c(Names[!isDF],DFNames)[duplicated(c(Names[!isDF],DFNames))])
	}
	for(i in seq_along(Names)){
		Sens <- ArgList[[i]]
		if(isNull[Names[i]]){
			Sname <- paste0("Sensor", Names[i])
		} else {
			Sname <- Names[i]
		}
		if(is.data.frame(Sens)){
			if(inherits(Sens, "Sensors") & length(Sens) == 8){
				Lnams <- names(Sens) <- c("name", "id", "node", "x", "y", "z", "d", "n")
			} else {
				Lnams <- names(Sens)
			}
			Sens <- as.data.frame(Sens)
			if(!all(c("x","y","z") %in% Lnams))stop(paste0("Sensor ",Names[i],": Please provide data.frame entries \"x\", \"y\" and \"z\"!\n"))
			if(is.character(Sens[1, 1])){
				Sens <- split(Sens, Sens[, 1])
			} else if(any(Sn <- (c("name", "sensor", "Sensor Name", "sensor_name") %in% Lnams))){
				Sens <- split(Sens, Sens[, c("name", "sensor", "Sensor Name", "sensor_name")[Sn]])
			} else {
				Sens <- setNames(list(Sens), Sname)
			}
		} else {
			if(any(is.na(Sens)))stop(paste0("Supplied list (argument ",Names[i],") contains NA values!\n"))
			Lnams <- names(Sens)
			if(!all(c("x","y","z") %in% Lnams))stop(paste0("Sensor ",Names[i],": Please provide list entries \"x\", \"y\" and \"z\"!\n"))
			Sens <- setNames(list(Sens), Sname)
		}

		for(snm in names(Sens)){
			dfin <- data.table(
				name = snm,
				x = as.numeric(Sens[[snm]][["x"]]),
				y = as.numeric(Sens[[snm]][["y"]]),
				z = as.numeric(Sens[[snm]][["z"]])
				)

			if("d" %in% Lnams){
				dfin[, d := as.numeric(Sens[[snm]][["d"]])]
			} else {
				# default = 1 m
				dfin[, d := 1.0]
			}

			if("id" %in% Lnams){
				dfin[, id := as.character(Sens[[snm]][["id"]])]
			} else {
				dfin[, id := "1"]
			}

			if("node" %in% Lnams){
				dfin[, node_orig := as.integer(Sens[[snm]][["node"]])]
			} else {
				dfin[, node_orig := seq_len(.N), by = id]
			}

			dfin[, z := round(z, 3)]

			dfin <- unique(dfin, by = c("name", "id", "x", "y", "z"))[, {
				data.table(node = order(node_orig), x, y, z, d)[order(node)]
			}, by = .(name, id)]

			dfout <- dfin[, {
				if(.N == 1L){
					.(node = 1L, x = x, y = y, z = z, d = 0.0, n = 1L)
				} else {
					N <- .N - 1
					dnew <- d * 0
					new <- dnew + 1
					for(j in seq_len(N)){
						l <- sqrt((x[j] - x[j+1]) ^ 2 + (y[j] - y[j+1]) ^ 2 + (z[j] - z[j+1]) ^ 2)
						if(d[j] > l){
							n1 <- 1
						} else {
							n1 <- l / d[j]
						}
						n <- round(n1 + 1, 0)
						dnew[j] <- round(l / (n - 1), 3)
						new[j] <- (n - 1)
					}
					.(node = node, x = x, y = y, z = z, d = dnew, n = new)
				}
			}, by = .(name, id)][, orig_i := i]
			Out <- rbind(Out, dfout)
		}
	}

	### check combined sensors
	if(Dups){
		for(dup in dupNames){
			Out[name == dup, id := {
				paste(id, frank(orig_i, ties.method = "dense"), sep = "_")
			}]
		}
	}
	Out[, orig_i := NULL]

	setnames(Out, c("Sensor Name", "Sensor ID", "Node", "x-Coord (m)", 
		"y-Coord (m)", "Sensor Height (m)", "Distance between Point-Sensors (m)", 
		"Number of Point-Sensors"))

	structure(as.data.frame(Out), class = c("Sensors", "data.frame"), Version = "4.2+")
}

# subset Sensors by name
'[.Sensors' <- function(x, i, j, ...) {
    if (!missing(i) && is.character(i)) {
        i <- which(x[[1]] %in% i)
    }
    out <- `[.data.frame`(x, i, j, ...)
    if (ncol(out) != 8 && inherits(out, 'data.frame')) {
        class(out) <- 'data.frame'
    }
    out
}
'[<-.Sensors' <- function(x, i, j, value) {
    if (!missing(i) && is.character(i)) {
        i <- which(x[[1]] %in% i)
    }
    out <- `[<-.data.frame`(x, i, j, value)
    if (ncol(out) != 8 && inherits(out, 'data.frame')) {
        class(out) <- 'data.frame'
    }
    out
}
