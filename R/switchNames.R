switchNames <- function(DT,simple=TRUE){
	sNames <- c(
		"rn","Sensor","Source","SensorHeight","SourceArea","CE","CE_se",
		"CE_lo","CE_hi","uCE","uCE_se","uCE_lo","uCE_hi","vCE","vCE_se","vCE_lo",
		"vCE_hi","wCE","wCE_se","wCE_lo","wCE_hi","N_TD","TD_Time_avg",
		"TD_Time_max","Max_Dist","UCE","Subset_seed","z_sWu","Ustar","L","Zo","sUu",
		"sVu","sWu","WD","d","kv","A","MaxFetch","N0","alpha","bw","C0","Sensor_Swustar"
		) 
	cNames <- c(
		"Original Interval Row","Sensor","Source",
		"Sensor Height [m above d]","Source Area [m]",
		"C/E","C/E SE","C/E Lower-95%CI","C/E Upper-95%CI",
		"u'C'/E","u'C'/E SE","u'C'/E Lower-95%CI","u'C'/E Upper-95%CI",
		"v'C'/E","v'C'/E SE","v'C'/E Lower-95%CI","v'C'/E Upper-95%CI",
		"w'C'/E","w'C'/E SE","w'C'/E Lower-95%CI","w'C'/E Upper-95%CI",
		"Number of TD inside","TD Time avg [sec]","TD Time max [sec]",
		"Maximum Distance [m]","uC/E","Catalog Subset Seed","SigmaW/Ustar Height [m]",
		"Ustar [m/s]","L [m]","Zo [m]","SigmaU/Ustar [-]","SigmaV/Ustar [-]",
		"SigmaW/Ustar [-]","WD [deg N]","d [m]","kv","A","MaxFetch","N0","alpha",
		"bw","C0","SigmaW/Ustar @ Sensor Height"
		)
	if(simple){
		from <- cNames
		to <- sNames
	} else {
		from <- sNames
		to <- cNames
	}
	UseMethod("switchNames",DT)
}
switchNames.data.table <- function(DT,simple=TRUE){
    from <- get('from', sys.frame(-1))
    to <- get('to', sys.frame(-1))
	cn <- names(DT)
	ex <- from %chin% cn
	setnames(DT,from[ex],to[ex])
	invisible(DT)
}
switchNames.data.frame <- function(DT,simple=TRUE){
    from <- get('from', sys.frame(-1))
    to <- get('to', sys.frame(-1))
	cn <- names(DT)
	ex <- chmatch(from,cn,nomatch=0)
	names(DT)[ex] <- to[as.logical(ex)]
	return(DT)
}


