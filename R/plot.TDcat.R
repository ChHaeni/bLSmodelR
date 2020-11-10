plot.TDcat <- function(x,asp=1,panel.first={grid();abline(h=0,col="darkgrey");abline(v=0,col="darkgrey")},...){
	invisible(x[,plot(x,y,asp=asp,panel.first=panel.first,...)])
}
