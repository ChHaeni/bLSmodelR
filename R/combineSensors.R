combineSensors <- function(res, comb_list = NULL, conc_column = NULL, add = FALSE){

  # convert old versions 
  sres <- as.character(substitute(res))
  res <- copy(res)
  setDT(res)
  if(is.null(attr(res, "Version"))){
    warning(paste0("Object '", sres[min(length(sres), 2)], "' has not yet been converted to version 4.2+"))
    convert(res)
  }

  if(is.null(comb_list)){
    uSensors <- res[, unique(Sensor)]
    comb_list <- setNames(list(uSensors), "all_sensors_combined")
  }

  # check duplicates
  if(any(duplicated(c(names(comb_list), res[, unique(Sensor)])))) stop("please provide new unique names for combined sensors")
  
  ModelInput <- attr(res, "ModelInput")

  outlist <- setNames(vector("list", length(comb_list)), names(comb_list))
  for(sensor_to in names(comb_list)){
    # sensor_to <- names(comb_list)[1]
    sensors <- comb_list[[sensor_to]]
    # get new sensor heights
    Sheight <- paste(sprintf("%1.2f", 
      range(as.numeric(unlist(strsplit(res[Sensor %in% sensors, unique(SensorHeight)], " to "))))
      ), collapse = " to ")
    # get Sensors
    Sensors <- ModelInput$Sensors[ModelInput$Sensors[, "Sensor Name"] %in% sensors, ]
    # check for point sensors:
    sn <- as.data.table(Sensors[, c("Sensor Name", "Sensor ID")])
    setnames(sn, c("name", "id"))
    nos <- sn[, .N, by = .(name, id)]
    if(nos[, any(N == 1)]){
      stop(paste0(
        "\nFollowing Sensor / ID are point sensors and cannot be combined: \n\t",
        paste(nos[N == 1, paste(name, id, sep = " / ")], collapse = "\n\t")
        )
      )
    }
    # get path lengths
    path_lengths <- getPathLengths(Sensors)
    outlist[[sensor_to]] <- res[Sensor %in% sensors, {
      if(.N == length(sensors)){
        out <- data.table(
          # Sensor
          Sensor = sensor_to,
          # SensorHeight
          SensorHeight = Sheight,
          # SourceArea
          SourceArea = SourceArea[1],
          # CE
          CE = avgCE_sensors(CE, path_lengths),
          # CE_se
          CE_se = avgCE_sensors_se(CE_se, path_lengths),
          CE_lo = NA_real_,
          CE_hi = NA_real_, 
          # uCE
          uCE = avgCE_sensors(uCE, path_lengths),
          # uCE_se
          uCE_se = avgCE_sensors_se(uCE_se, path_lengths),
          uCE_lo = NA_real_,
          uCE_hi = NA_real_, 
          # vCE
          vCE = avgCE_sensors(vCE, path_lengths),
          # vCE_se
          vCE_se = avgCE_sensors_se(vCE_se, path_lengths),
          vCE_lo = NA_real_, 
          vCE_hi = NA_real_,
          # wCE
          wCE = avgCE_sensors(wCE, path_lengths),
          # wCE_se
          wCE_se = avgCE_sensors_se(wCE_se, path_lengths),
          wCE_lo = NA_real_, 
          wCE_hi = NA_real_,
          # N_TD
          N_TD = sum(N_TD),
          # TD_Time_avg
          TD_Time_avg = avgCE_sources(TD_Time_avg, N_TD / sum(N_TD)),
          # TD_Time_max
          TD_Time_max = max(TD_Time_max),
          # Max_Dist
          Max_Dist = max(Max_Dist),
          # UCE
          UCE = avgCE_sensors(UCE, path_lengths)
          )
        # all lo und hi nachrechnen...
        N0m <- mean(N0)
        qtlo <- qt(0.025, N0m - 1)
        qthi <- qt(0.975, N0m - 1)
        out[, ":="(
          CE_lo = CE + qtlo * CE_se,
          CE_hi = CE + qthi * CE_se,
          uCE_lo = uCE + qtlo * uCE_se,
          uCE_hi = uCE + qthi * uCE_se,
          vCE_lo = vCE + qtlo * vCE_se,
          vCE_hi = vCE + qthi * vCE_se,
          wCE_lo = wCE + qtlo * wCE_se,
          wCE_hi = wCE + qthi * wCE_se                 
          )]
        # add concentration column if any
        if(!is.null(conc_column)){
          out[, (conc_column) := avgCE_sensors(get(conc_column), path_lengths)]
        }
        # cbind all missing columns!
        cbind(
          out,
          .SD[1, names(.SD) %w/o% names(out), with = FALSE]
          )
      }
    }, by = .(rn, Source)]
  }

  # adjust attributes part 1
  CalcSteps <- attr(res, "CalcSteps")
  ModelInput <- attr(res, "ModelInput")
  Catalogs <- attr(res, "Catalogs")
  ModelInput$Sensors <- NULL
  CalcSteps_out <- NULL
  Catalogs_out <- NULL

  # add sensors to Sensor
  for(sensors_to in names(comb_list)){
    # sensors_to <- names(comb_list)[1]
    sensors <- comb_list[[sensors_to]]
    # get all point sensors
    Sensors <- attr(res, "ModelInput")$Sensors[attr(res, "ModelInput")$Sensors[, "Sensor Name"] %in% sensors, ]
    # old sensors
    Calc.Sensors <- procSensors(Sensors)$Calc.Sensors
    point_sensors <- Calc.Sensors[, "Point Sensor Name"]
    # new sensors
    NewSensors <- Sensors
    # IDs & Names
    NewSensors[, "Sensor ID"] <- paste(Sensors[, "Sensor Name"], Sensors[, "Sensor ID"], sep = "_")
    NewSensors[, "Sensor Name"] <- sensors_to
    NewCalc.Sensors <- procSensors(NewSensors)$Calc.Sensors
    new_point_sensors <- NewCalc.Sensors[, "Point Sensor Name"]
    names(new_point_sensors) <- point_sensors
    # get rows
    Rows <- CalcSteps[, .(take = all(sensors %in% Sensor)), by = .(rn, Source)]
    # CalcSteps
    CalcSteps_temp <- merge(CalcSteps, Rows[(take), .(rn, Source)], by = c("rn", "Source"))
    # Catalogs
    Catalogs_temp <- merge(Catalogs, Rows[(take), .(rn)], by = "rn")[PointSensor %chin% point_sensors]
    # replace CalcSteps
    CalcSteps_temp[, Sensor := sensors_to]
    CalcSteps_temp[, New.Sensor := {
      paste(
        new_point_sensors[match(point_sensors, 
          unlist(strsplit(.BY$Old_Sensor, ",")), nomatch = 0)],
        collapse = ",")
    }, by = .(Old_Sensor = Calc.Sensor)][, ":="(
      Calc.Sensor = New.Sensor,
      New.Sensor = NULL
      )]
    # replace Catalogs
    Catalogs_temp[, Sensor := sensors_to]
    Catalogs_temp[, New.Sensor := {
      new_point_sensors[.BY$Old_Sensor]
    }, by = .(Old_Sensor = PointSensor)][, ":="(
      PointSensor = New.Sensor,
      New.Sensor = NULL
      )]
    # add
    ModelInput$Sensors <- rbind(
      ModelInput$Sensors,
      NewSensors
      )
    CalcSteps_out <- rbind(
      CalcSteps_out,
      CalcSteps_temp[Calc.Sensor != ""]
      )
    Catalogs_out <- rbind(
      Catalogs_out,
      Catalogs_temp
      )
  }

  if(add){
    out <- rbind(
      res,
      rbindlist(outlist)
      )
    # add
    ModelInput$Sensors <- rbind(
      attr(res, "ModelInput")$Sensors,
      ModelInput$Sensors
      )
    CalcSteps_out <- rbind(
      CalcSteps_out,
      attr(res, "CalcSteps")
      )
    Catalogs_out <- rbind(
      Catalogs_out,
      attr(res, "Catalogs")
      )
  } else {
    out <- rbindlist(outlist)   
  }
  setorder(CalcSteps_out, rn, Sensor, Source)

  setkey(CalcSteps_out, rn, Sensor)
  setkey(Catalogs_out, rn, Sensor, PointSensor)

  # adjust attributes part 2
  for(att in names(attributes(res)) %w/o% c("row.names", 
    ".internal.selfref", "names", "index", "Catalogs", "CalcSteps", "ModelInput")){
    setattr(out, att, attr(res, att))
  }
  setattr(out, "CalcSteps", CalcSteps_out)
  setattr(out, "Catalogs", Catalogs_out)
  setattr(out, "ModelInput", ModelInput)
  setattr(out, "combined_sensors_list", comb_list)
  
  setorder(out, rn, Sensor, Source)
  setindex(out, Source)

  out
}

avgCE_sensors <- function(ce, paths){
  sum(paths * ce) / sum(paths)
}
avgCE_sensors_se <- function(ce_se, paths){
  sqrt(sum(paths ^ 2 * ce_se ^ 2)) / sum(paths)
}

getPathLengths <- function(x){
  snsrs <- as.data.table(x[, c("Sensor Name", "Sensor ID", "Node", "x-Coord (m)", "y-Coord (m)", "Sensor Height (m)")])
  setnames(snsrs, c("name", "id", "node", "x", "y", "z"))
  snsrs[,{
    sum(sqrt(diff(x) ^ 2 + diff(y) ^ 2 + diff(z) ^ 2))
  }, by = .(name, id)][,{
    sum(V1)
  }, by = name][, setNames(V1, name)]
}
