coreModelWrapper <- function(index,u,v,w,sncrun){
	coreModel(u[index],v[index],w[index],sncrun[,SensorHeight],sncrun[,Ustar],sncrun[,L],sncrun[,Zo],sncrun[,bw],sncrun[,sUu],sncrun[,sVu],sncrun[,kv],sncrun[,C0],sncrun[,alpha],sncrun[,MaxFetch])
}
