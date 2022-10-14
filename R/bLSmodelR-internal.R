.onLoad <- function(libname, pkgname) {
	 invisible(
         reg.finalizer(
        e = parent.env(environment()),
        f = function(env){
            eval(env, {
                bLSmodelR::cleanTemporary()
            })
        },
        onexit = TRUE))
}
.onAttach <- function(libname, pkgname){
    packageStartupMessage("\n#################################")
    packageStartupMessage(paste0(" This is bLSmodelR version ",packageVersion("bLSmodelR")))
    packageStartupMessage(paste0(" last updated ",packageDescription("bLSmodelR")[["Date"]]))
    packageStartupMessage("#################################\n")
    # set some options
    options(
        # parent directory of job dir
        bls.slurm.jobdir = getOption('bls.slurm.jobdir', file.path(Sys.getenv('HOME'), '.slurm')),
        # exclude partitions?
        bls.slurm.exclude.partition = getOption('bls.slurm.exclude.partition', '')
    )
}
.inifu <- function(n,theta){
	# http://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variable
	Id <- diag(n)
	x1    <- rnorm(n)        # fixed given data
	x1 <- (x1 - mean(x1))#/sd(x1)
	x2    <- rnorm(n)      # new random data
	Xctr   <- cbind(x1, x2 - mean(x2))
	Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
	P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
	x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
	Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
	Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
	x3    <- rnorm(n)      # new random data
	Xctr     <- cbind(Y, x3 - mean(x3))         # matrix
	Q    <- qr.Q(qr(Xctr[ , 1:2, drop=FALSE]))      # QR-decomposition, just matrix Q
	P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
	x3o  <- (Id-P) %*% Xctr[ , 3]                 # x2ctr made orthogonal to x1ctr
	Xc2  <- cbind(Xctr[ , 1], x3o)                # bind to matrix
	Y2    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
	x2 <- Y[ , 2]    # final new vector
	x3 <- Y2[ , 2] + (1 / tan(theta)) * Y[ , 1]    # final new vector
	return(cbind(x1/sd(x1),x2/sd(x2),x3/sd(x3)))
}
.MatchWrapper <- function(Int_Ext,Cat_list,Tol_){
	# check matching
	Int_Ext[,{
		# browser()
		Cat_list[
			abs(Cat_kv - .BY$kv) < 1E-2 &
			abs(Cat_A - .BY$A) < 1E-2 &
			abs(Cat_alpha - .BY$alpha) < 1E-3 &
			Cat_MaxFetch >= .BY$MaxFetch &
			Cat_ZSens <= SensorHeight_Upper &
			Cat_ZSens >= SensorHeight_Lower &
			sign(Cat_L) == sign(.BY$L) &
			abs(Cat_L) <= L_Upper &
			abs(Cat_L) >= L_Lower &
			Cat_Zo <= Zo_Upper &
			Cat_Zo >= Zo_Lower &
			Cat_Su_Ustar <= sUu_Upper &
			Cat_Su_Ustar >= sUu_Lower &
			Cat_Sv_Ustar <= sVu_Upper &
			Cat_Sv_Ustar >= sVu_Lower &
			Cat_Sensor_Swustar <= sWu_Upper &
			Cat_Sensor_Swustar >= sWu_Lower
		,.SD]
		},by = .(rn,Cat.Name,SensorHeight,L,Zo,sUu,sVu,Sensor_Swustar,MaxFetch,Sensor,N0,alpha,A,kv)][
			,":="(
			devZSens = abs(Cat_ZSens/SensorHeight - 1)/Tol_[1],
			devL = abs(Cat_L/L - 1)/Tol_[2],
			devZo = abs(Cat_Zo/Zo - 1)/Tol_[3],
			devsUu = abs(Cat_Su_Ustar/sUu - 1)/Tol_[4],
			devsVu = abs(Cat_Sv_Ustar/sVu - 1)/Tol_[5],
			devSensor_Swustar = abs(Cat_Sensor_Swustar/Sensor_Swustar - 1)/Tol_[6],
			devN0 = Cat_N0 - N0
			)
		][
			devZSens <= 1.0000001 &
			devL <= 1.0000001 &
			devZo <= 1.0000001 &
			devsUu <= 1.0000001 &
			devsVu <= 1.0000001 &
			devSensor_Swustar <= 1.0000001,sumDev:=devZSens+devL+devZo+devsUu+devsVu+devSensor_Swustar]

}
.CrossMatchWrapper <- function(Int_Ext,Cat_list,Tol_){

			# check matching heights etc.
			CheckPart1 <- Int_Ext[, .(rn, SensorHeight,Sensor_Swustar,MaxFetch,
				sWu_Upper,
				sWu_Lower,
				SensorHeight_Upper,
				SensorHeight_Lower
				)][,{
				Cat_list[
					Cat_MaxFetch >= MaxFetch[1] &
					Cat_Sensor_Swustar <= sWu_Upper[1] &
					Cat_Sensor_Swustar >= sWu_Lower[1] &
					Cat_ZSens <= SensorHeight_Upper[1] &
					Cat_ZSens >= SensorHeight_Lower[1],{
					.(
						SensorHeight = SensorHeight[1], 
						Sensor_Swustar = Sensor_Swustar[1], 
						MaxFetch = MaxFetch[1], 
						Zeile = Zeile)
				}]		
			}, keyby = .(ssm = paste(SensorHeight, Sensor_Swustar, MaxFetch, sep = "/"))]


			# check matching MOST
			Int_Ext[,.(rn,Cat.Name,SensorHeight,L,Zo,sUu,sVu,Sensor_Swustar,MaxFetch,Sensor,N0,
				sUu_Upper,
				sUu_Lower,
				sVu_Upper,
				sVu_Lower,
				Zo_Upper,
				Zo_Lower,
				L_Upper,
				L_Lower
				)][,{

				ind <- paste(SensorHeight, Sensor_Swustar, MaxFetch, sep = "/")
				Sub <- CheckPart1[.(ind), .(Zeile, ssm)]
				ind2 <- Sub[,unique(Zeile)]
				out <- Cat_list[match(ind2,Zeile)][
					Cat_Su_Ustar <= sUu_Upper[1] &
					Cat_Su_Ustar >= sUu_Lower[1] &
					Cat_Sv_Ustar <= sVu_Upper[1] &
					Cat_Sv_Ustar >= sVu_Lower[1] &
					Cat_Zo <= Zo_Upper[1] &
					Cat_Zo >= Zo_Lower[1] &
					sign(Cat_L) == sign(L[1]) &
					abs(Cat_L) <= L_Upper[1] &
					abs(Cat_L) >= L_Lower[1],
						merge(.SD, Sub, by = "Zeile")]

				c(
					.SD[out[,match(ssm, ind)],],
					out
					)
			}, by = rn][
				,":="(
				devZSens = abs(Cat_ZSens/SensorHeight - 1)/Tol_[1],
				devL = abs(Cat_L/L - 1)/Tol_[2],
				devZo = abs(Cat_Zo/Zo - 1)/Tol_[3],
				devsUu = abs(Cat_Su_Ustar/sUu - 1)/Tol_[4],
				devsVu = abs(Cat_Sv_Ustar/sVu - 1)/Tol_[5],
				devSensor_Swustar = abs(Cat_Sensor_Swustar/Sensor_Swustar - 1)/Tol_[6],
				devN0 = Cat_N0 - N0
				)
			][
				devZSens <= 1.0000001][
				devL <= 1.0000001][
				devZo <= 1.0000001][
				devsUu <= 1.0000001][
				devsVu <= 1.0000001][
				devSensor_Swustar <= 1.0000001,sumDev:=devZSens+devL+devZo+devsUu+devsVu+devSensor_Swustar]

}

.MaxFetchWrapper <- function(Int_Ext, p_Sens, Input_List){
	Int_Ext[MaxFetch < 0,
		MaxFetch := {
			nmSens <- unlist(strsplit(Sensor,split=","))
			nmSous <- unlist(strsplit(Source,split=","))
			Sens <- structure(p_Sens$Calc.Sensors[p_Sens$Calc.Sensors[, "Point Sensor Name"] %in% nmSens,
				c("x-Coord (m)", "y-Coord (m)")],class="data.frame")
			Sous <- structure(Input_List[["Sources"]][Input_List[["Sources"]][,1] %in% nmSous,2:3],class="data.frame")
			out <- numeric(length(WD))
			for(i in seq_along(WD)){
				Sensrot <- rotate(Sens,-WD[i])
				Sourot <- rotate(Sous,-WD[i])
				out[i] <- max(Sensrot[,1]) - min(Sourot[,1])  
			}
			out <- ceiling(pmax(out,0)) - MaxFetch
			out
		}
	,by=.(Sensor,Source)]		
}
