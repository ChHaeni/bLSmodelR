createCatName <- function(dt,Tstring=""){
	dt[,
	{ZoString <- sub("[.]","",sub("e[+]","Ep",sub("e[-]","Em",sprintf("%1.3e",Zo))))
	LString <- gsub("[.-]","",sub("e[+]","Ep",sub("e[-]","Em",sprintf("%1.3e",L))))
	sprintf("Cat_Zm%1.0f_L%s%s_Zo%s_bw%1.0f_%s",
		1000*SensorHeight,
		ifelse(L<0,"m","p"),
		LString,
		ZoString,
		1000*bw,
		Tstring)
	}]
}
