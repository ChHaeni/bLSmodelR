

procSensors <- function(x){
	nms <- names(x)
	x <- as.data.table(x)
	setnames(x, c("name", "id", "node", "x", "y", "z", "d", "n"))
	x <- x[,{
		data.table(name, new_id = order(id), node, x, y, z, d, n)[order(new_id)]
	}, by = .(name, id)]
	## extend paths to individual point sensors
	Snsrs <- x[,{
		N <- .N - 1
		out <- data.table(node = 0L, x = x[1], y = y[1], z = z[1])
		for(i in seq_len(N)){
			sub_ind <- (as.numeric(i != 1) + 1):(n[i] + 1)
			xout <- seq(x[i], x[i + 1], length.out = n[i] + 1)
			yout <- seq(y[i], y[i + 1], length.out = n[i] + 1)
			zout <- round(seq(z[i], z[i + 1], length.out = n[i] + 1), 3)
			out <- rbind(out, data.table(node = 0L, x = xout[-1], y = yout[-1], z = zout[-1]))
		}
		out[, node := seq_len(.N)]
	}, by = .(name, id)]
	Snsrs[, point_sensor := paste(name, id, node, sep = ".")]

	### get heights
	heights <- Snsrs[, unique(z)]
	hts <- lapply(heights,function(x,y){y[y[,2]==x,1]},y=Snsrs[, data.frame(point_sensor,z, stringsAsFactors = FALSE)])
	names(hts) <- heights

	# sort by height:
	sindex <- order(as.numeric(names(hts)),decreasing=TRUE)
	hts <- hts[sindex]

	# Line/path Sensor Names
	Lnames <- Snsrs[, {
		tNames <- table(name)
		names(tNames)[tNames > 1]
	}]

	# set names
	setnames(Snsrs, c(nms[1:6], "Point Sensor Name"))
	setDF(Snsrs)

	# calc PS_list
	PS_list <- tapply(Snsrs[, "Point Sensor Name"], Snsrs[, "Sensor Name"], function(x){
		paste(x, collapse = ",")
	}, simplify = FALSE)

	list(
		heights = hts, 
		Calc.Sensors = Snsrs, 
		LineSensors = Lnames,
		PS_list = PS_list
		)
}
