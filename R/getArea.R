getArea <- function(x,individual.polys=FALSE){
	if(inherits(x,"Sources") && ncol(x) == 4){
		Srcs <- unique(x[,1])
		if(individual.polys){
			out <- vector("list",length=length(Srcs))
			names(out) <- Srcs
			for(i in Srcs){
				y <- x[x[[1]] %in% i,]
				polys <- unique(y[,4])
				out[[i]] <- numeric(length(polys))
				names(out[[i]]) <- polys
				for(j in polys){
					z <- y[y[[4]] %in% j,2:3]
					ind <- c(nrow(z), seq(nrow(z) - 1))
					out[[i]][j] <- abs(sum(z[,1]*z[ind,2] - z[ind,1]*z[,2]))/2
				}
			}			
		} else {
			out <- numeric(length(Srcs))
			names(out) <- Srcs
			for(i in Srcs){
				y <- x[x[[1]] %in% i,]
				polys <- unique(y[,4])
				for(j in polys){
					z <- y[y[[4]] %in% j,2:3]
					ind <- c(nrow(z), seq(nrow(z) - 1))
					out[i] <- out[i] + abs(sum(z[,1]*z[ind,2] - z[ind,1]*z[,2]))/2
				}
			}			
		}
		return(out)
	} else {
	    ind <- c(nrow(x), seq(nrow(x) - 1))
	    return(abs(sum(x[,1]*x[ind,2] - x[ind,1]*x[,2]))/2)
	}
}
