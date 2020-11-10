print.TDcat <- function(x,...){
	x <- copy(x)
	setattr(x,"class",c("data.table","data.frame"))
	if(!is.null(attr(x,"N0"))){
		cat("\n")
		cat(attr(x,"header"))
		cat("\nTouchdown Catalog:\n")
		if(nrow(x)>0){
			ucdim1 <- unlist(strsplit(sprintf("%d",attr(x,"N0")),""))
			add <- c(rep(c("","","'"),floor((length(ucdim1)-1)/3)),rep("",3))
			pdim1 <- paste(ucdim1,paste(add[(1:length(ucdim1))+(length(add)-length(ucdim1))],sep=""),sep="",collapse="")
			cat(paste0("\n",pdim1," trajectories:\n"))
			ucdim1 <- unlist(strsplit(sprintf("%d",nrow(x)),""))
			add <- c(rep(c("","","'"),floor((length(ucdim1)-1)/3)),rep("",3))
			pdim1 <- paste(ucdim1,paste(add[(1:length(ucdim1))+(length(add)-length(ucdim1))],sep=""),sep="",collapse="")
			cat(paste0("\n",pdim1," touchdowns:\n"))
			print(x)
		} else {
			cat("\nempty TD Catalog.\n\n")
		}
		cat("\ninitialized velocities:\n")
		if(!is.null(attr(x,"uvw0"))){
			if(attr(x,"N0")>6){
				hx <- head(attr(x,"uvw0"),3)
				tx <- tail(attr(x,"uvw0"),3)
				px <- rbind(hx,rep(NA,3),tx)
				rownames(px) <- c(paste("[",1:3,", ]",sep=""),"***",rownames(tx))
				print(px,na.print="***",...)			
			} else {
				print(attr(x,"uvw0"))
			}		
		} else {
			cat("\ntriggered... (run \"initializeCatalog\" to finish)\n")
		}
		cat("\n")
	} else {
		print(x)
	}
}
