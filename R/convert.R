
convert <- function(x){
  if(is.null(attr(x, "Version"))){
    if(inherits(x, "bLSresult")){
      isdt <- is.data.table(x)
      setDT(x)
      CalcSteps <- attr(x, "CalcSteps")
      Catalogs <- attr(x, "Catalogs")
      ModelInput <- attr(x, "ModelInput")
      # add seed
      CalcSteps[x[, .(rn, Sensor, Subset_seed)], seed := Subset_seed, on = c("rn", "Sensor")]
      Catalogs[x[, .(rn, Sensor, Subset_seed)], seed := Subset_seed, on = c("rn", "Sensor")]
      # change Sensors
      ModelInput$Sensors <- convert(ModelInput$Sensors)
      # change SensorHeight
      Sensors <- as.data.table(ModelInput$Sensors)
      setnames(Sensors, c("Sensor Name", "Sensor Height (m)"), c("name", "z"))
      Sht <- Sensors[, {
        if(.N > 1){
          paste(sprintf("%1.2f", range(z)), collapse = " to ")
        } else {
          sprintf("%1.3f", z)
        }
      }, by = name][, setNames(V1, name)]
      x[, SensorHeight := as.character(SensorHeight)]
      x[, SensorHeight := Sht[Sensor]]

      # add seed
      CSens <- procSensors(ModelInput$Sensors)$Calc.Sensors
      old <- CSens[, "Sensor Name"]
      numS <- table(CSens[, "Sensor Name"])
      ind <- old %in% names(numS)[numS > 1]
      old[ind] <- paste(old[ind], CSens[, "Node"], sep = ".")
      new <- setNames(CSens[, "Point Sensor Name"], old)
      CalcSteps[, Calc.Sensor := {
        sens <- unlist(strsplit(.BY[[1]], split = ","))
        paste(new[sens], collapse = ",")
      }, by = Calc.Sensor]
      Catalogs[, PointSensor := new[PointSensor]]

      setattr(x, "CalcSteps", CalcSteps)
      setattr(x, "Catalogs", Catalogs)
      setattr(x, "ModelInput", ModelInput)
      # remove columns
      x[, c("Sensor_Swustar", "Calc.ZSens", "Calc.Ustar", "Calc.L", "Calc.Zo",
        "Calc.Su_Ustar", "Calc.Sv_Ustar", "Calc.bw", "Calc.C0", "Calc.kv", 
        "Calc.A", "Calc.alpha", "Calc.MaxFetch", "Calc.Sensor_Swustar", 
        "Calc.N0", "Subset_seed") := NULL]
      setattr(x, "Version", "4.2+")
      if(!isdt){
        setDF(x)
      }
      return(invisible(x))
    } else if(inherits(x, "Sensors")){
      # change Sensors
      Sensors <- x[, c(1:4, 7)]
      Sensors[x[,5] != "", 1] <- x[x[,5] != "", 5]
      names(Sensors) <- c("name", "x", "y", "z", "d")
      x <- genSensors(Sensors)
      attr(x, "Version") <- "4.2+"    
      return(x)
    } else if(inherits(x, "InputList")){
      # change Sensors
      x$Sensors <- convert(x$Sensors)  
      return(x)
    }
  } else {
    invisible(x)
  } 
}
