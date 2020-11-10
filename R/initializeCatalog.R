initializeCatalog <- function(N0,Ustar,Suu,Svu,bw,SensorHeight,L,Zo,A,alpha,MaxFetch,kv,Catalog=NULL){
	env <- globalenv()
    oseed <- env$.Random.seed

	if(is.null(Catalog)){
		if(is.data.frame(N0)){
			DT <- as.data.table(N0)
			if(nrow(DT)!=1)stop("Please provide only 1 single row!")
			N0 <- DT[,as.numeric(N0)]
			Ustar <- DT[,as.numeric(Ustar)]
			Suu <- DT[,as.numeric(sUu)]
			Svu <- DT[,as.numeric(sVu)]
			bw <- DT[,as.numeric(bw)]
			SensorHeight <- DT[,as.numeric(SensorHeight)]
			L <- DT[,as.numeric(L)]
			Zo <- DT[,as.numeric(Zo)]
			A <- DT[,as.numeric(A)]
			alpha <- DT[,as.numeric(alpha)]
			MaxFetch <- DT[,as.numeric(MaxFetch)]
			kv <- DT[,as.numeric(kv)]	
		}
		U0 <- calcU(Ustar, Zo, L, SensorHeight, kv)
		Catalog <- data.table(Traj_ID=integer(0),Time=numeric(0),x=numeric(0),y=numeric(0),wTD=numeric(0))
		setattr(Catalog,"Ustar",Ustar) 
		setattr(Catalog,"U0",U0)
		setattr(Catalog,"N0",N0)
		setattr(Catalog,"UVWseed",sample.int(1E9,1))
		setattr(Catalog,"class",c("TDcat","data.table","data.frame"))
		txt <- sprintf("N0 = %i\nZSens = %1.3f\nUstar = %1.7f\nL = %1.1f\nZo = %1.3e\nSu_Ustar = %1.3f\nSv_Ustar = %1.3f\nbw = %1.3f\nC0 = %1.3f\nkv = %1.2f\nA = %1.2f\nalpha = %1.3f\nMaxFetch = %1.0f\n",N0,SensorHeight,Ustar,L,Zo,Suu,Svu,bw,calcC0(bw,kv,A),kv,A,alpha,MaxFetch)
		TDhead <- sprintf(paste0("TD Catalog generated on %s\n",txt),format(Sys.time()))
		class(TDhead) <- c("TDhead","character")
		setattr(Catalog,"header",TDhead)
		# initializeUVW:
	    on.exit({if (is.null(oseed)) {
				rm(list = ".Random.seed", envir = env)
			} else {
				assign(".Random.seed", value = oseed, envir = env)
			}
		})			
		set.seed(attr(Catalog,"UVWseed"),kind="L'Ecuyer-CMRG")
		Su <- Suu*Ustar
		Sv <- Svu*Ustar
		Sw <- calcsigmaW(Ustar,SensorHeight/L,bw)
		theta  <- acos(-Ustar/(Suu*Sw))
		N2 <- max(floor(N0/100),1)
		n <- rep(100,N2)
		n[1] <- n[1] + N0%%100
		out <- matrix(NA,nrow=N0,ncol=3)
		dummy <- 0
		for(i in 1:N2){
			dummy <- max(dummy) + 1:n[i]
			out[dummy,] <- .inifu(n[i],theta)
		}
		out <- t(t(out)/apply(out,2,sd))		
		setattr(Catalog,"uvw0",cbind(u0=out[,1]*Su + U0,v0=out[,2]*Sv,w0=out[,3]*Sw))
		setattr(Catalog,"is.compact",FALSE)
	} else {
		if(is.data.frame(N0)){
			DT <- as.data.table(N0)
			if(nrow(DT)!=1)stop("Please provide only 1 single row!")
			N0 <- DT[,as.numeric(N0)]
			Ustar <- DT[,as.numeric(Ustar)]
		}
		hdr <- attr(Catalog,"header")
		SensorHeight <- as.numeric(gsub(".*ZSens = ([0-9]+[.][0-9]{3}).*","\\1",hdr))
		L <- as.numeric(gsub(".*L = ([-]?[0-9]+[.][0-9]{1}).*","\\1",hdr))
		Zo <- as.numeric(gsub(".*Zo = (.*).Su.*","\\1",hdr))
		A <- as.numeric(gsub(".*A = ([0-9]+[.][0-9]{2}).*","\\1",hdr))
		alpha <- as.numeric(gsub(".*alpha = ([0-9]+[.][0-9]{3}).*","\\1",hdr))
		MaxFetch <- as.numeric(gsub(".*MaxFetch = ([0-9]+).*","\\1",hdr))
		kv <- as.numeric(gsub(".*kv = ([0-9]+[.][0-9]{2}).*","\\1",hdr))		
		Suu <- sd(attr(Catalog,"uvw0")[,"u0"])/attr(Catalog,"Ustar")
		Svu <- sd(attr(Catalog,"uvw0")[,"v0"])/attr(Catalog,"Ustar")
		bw <- calcbw(sd(attr(Catalog,"uvw0")[,"w0"])/attr(Catalog,"Ustar"),SensorHeight/L)

		Cat_N0 <- attr(Catalog,"N0")
		Nresid <- N0 - Cat_N0
		if(Nresid>0){
			if(Nresid>99){
				# getUVW:
				# initializeUVW:
			    on.exit({if (is.null(oseed)) {
						rm(list = ".Random.seed", envir = env)
					} else {
						assign(".Random.seed", value = oseed, envir = env)
					}
				})			
				set.seed(attr(Catalog,"UVWseed"),kind="L'Ecuyer-CMRG")
				Su <- Suu*attr(Catalog,"Ustar")
				Sv <- Svu*attr(Catalog,"Ustar")
				Sw <- calcsigmaW(attr(Catalog,"Ustar"),SensorHeight/L,bw)
				theta  <- acos(-attr(Catalog,"Ustar")/(Suu*Sw))
				N2 <- max(floor(Cat_N0/100),1)
				n <- rep(100,N2)
				n[1] <- n[1] + Cat_N0%%100
				out <- matrix(NA,nrow=Cat_N0,ncol=3)
				dummy <- 0
				for(i in 1:N2){
					dummy <- max(dummy) + 1:n[i]
					out[dummy,] <- .inifu(n[i],theta)
				}
				out <- t(t(out)/apply(out,2,sd))	
				# add uvw:
				N2 <- max(floor(Nresid/100),1)
				n <- rep(100,N2)
				n[1] <- n[1] + Nresid%%100
				outadd <- matrix(NA,nrow=Nresid,ncol=3)
				dummy <- 0
				for(i in 1:N2){
					dummy <- max(dummy) + 1:n[i]
					outadd[dummy,] <- .inifu(n[i],theta)
				}
				outadd <- t(t(outadd)/apply(outadd,2,sd))			
				pars <- numeric(3)
				for(i in 1:3)pars[i] <- optimize(function(x)abs(1-sd(c(out[,i],outadd[,i]/x))),interval=c(0.5,1.5),tol=1E-10)$minimum
				out <- rbind(out,t(t(outadd)/pars))
				setattr(Catalog,"uvw0",cbind(u0=out[,1]*Su + attr(Catalog,"U0"),v0=out[,2]*Svu*attr(Catalog,"Ustar"),w0=out[,3]*Sw))
				setattr(Catalog,"N0",N0)
				txt <- sprintf("N0 = %i\nZSens = %1.3f\nUstar = %1.7f\nL = %1.1f\nZo = %1.3e\nSu_Ustar = %1.3f\nSv_Ustar = %1.3f\nbw = %1.3f\nC0 = %1.3f\nkv = %1.2f\nA = %1.2f\nalpha = %1.3f\nMaxFetch = %1.0f\n",N0,SensorHeight,attr(Catalog,"Ustar"),L,Zo,Suu,Svu,bw,calcC0(bw,kv,A),kv,A,alpha,MaxFetch)
				TDhead <- sprintf(paste0("TD Catalog generated on %s\n",txt),format(Sys.time()))
				class(TDhead) <- c("TDhead","character")
				setattr(Catalog,"header",TDhead)
				# check compact:
				if(attr(Catalog,"is.compact")){
					if(attr(Catalog,"is.int")){
	 					tdvalues <- 10^attr(Catalog,"tdvalues.round")
		            	Catalog[,":="(
		                    Time = as.numeric(Time/tdvalues[1])
		                    ,x = as.numeric(x/tdvalues[2])
		                    ,y = as.numeric(y/tdvalues[3])
		                    ,wTD = as.numeric(wTD/tdvalues[4])
		                )]
		                setattr(Catalog,"is.int",FALSE)
	 				}
	 				setattr(Catalog,"uvw.round",c("u0"=Inf,"v0"=Inf,"w0"=Inf))
				}
			} else {
				warning("N0: ",N0,"\nCatalog: ",Cat_N0,"\nIncrease N0 by at least 100 trajectories compared to the available Catalog trajectory number!\n-> N0 is reset to ",Cat_N0,"!\n")
			}
		} else {
			N0 <- Cat_N0
		}
		# Scale Ustar
		Uratio <- Ustar/attr(Catalog,"Ustar")
		if(Uratio!=1){
			Catalog[,wTD:=wTD*Uratio]
			Catalog[,Time:=Time/Uratio]
			setattr(Catalog,"uvw0",attr(Catalog,"uvw0")*Uratio)
			setattr(Catalog,"Ustar",Ustar)
			setattr(Catalog,"U0",attr(Catalog,"U0")*Uratio)
		}
	}
	
	return(Catalog)
}
