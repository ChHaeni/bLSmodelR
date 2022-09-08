combineSensors <- function(result, comb_list = NULL, conc_column = NULL, add = FALSE){

    # convert old versions 
    has_dep <- inherits(result, 'deposition')
    res <- copy(result)
    setDT(res)
    if(is.null(attr(res, "Version"))){
        sres <- as.character(substitute(result))
        warning(paste0("Object '", sres[min(length(sres), 2)], "' has not yet been converted to version 4.2+"))
        convert(res)
    }

    if(is.null(comb_list)){
        uSensors <- res[, unique(Sensor)]
        comb_list <- setNames(list(uSensors), "all_sensors_combined")
    }

    # check duplicates
    if(any(duplicated(c(names(comb_list), res[, unique(Sensor)])))) stop("please provide new unique names for combined sensors")

    # get attributes part 1
    ModelInput <- attr(res, "ModelInput")
    CalcSteps_out <- NULL
    Catalogs_out <- NULL

    # get unique sensors
    uSensors <- unique(unlist(comb_list))

    # check Sensors
    Sensors <- ModelInput$Sensors[ModelInput$Sensors[, "Sensor Name"] %in% uSensors, ]
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
    path_lengths <- getPathLengths(ModelInput$Sensors)[uSensors]
    comb_lengths <- lapply(comb_list, function(x) {
        tot <- sum(path_lengths[x])
        list(
            wts = path_lengths[x] / tot,
            wts_sqr = (path_lengths[x] / tot) ^ 2,
            tot = tot
            )
    })

    # subsets & attributes part 2
    res <- res[Sensor %in% uSensors]
    CalcSteps <- attr(res, 'CalcSteps')[Sensor %in% uSensors]
    Catalogs <- attr(res, "Catalogs")[Sensor %in% uSensors]
    # erase Sensors...
    ModelInput$Sensors <- NULL

    outlist <- setNames(vector("list", length(comb_list)), names(comb_list))
    for(sensor_to in names(comb_list)) {
        # get sensor names
        sensors_from <- comb_list[[sensor_to]]
        cat('Combining sensors', paste(sensors_from, collapse = ' + '), 'to', sensor_to, '\n')
        # get new sensor heights
        Sheight <- paste(sprintf("%1.2f", 
                range(as.numeric(unlist(strsplit(res[Sensor %in% sensors_from, unique(SensorHeight)], " to "))))
                ), collapse = " to ")
        # average values
        outlist[[sensor_to]] <- res[Sensor %chin% sensors_from, {
            if (.N == length(sensors_from)) {
                wts <- comb_lengths[[sensor_to]][['wts']][Sensor]
                wts_sqr <- comb_lengths[[sensor_to]][['wts_sqr']][Sensor]
                out <- data.table(
                    # Sensor
                    Sensor = sensor_to,
                    # SensorHeight
                    SensorHeight = Sheight,
                    # SourceArea
                    SourceArea = SourceArea[1],
                    # CE
                    CE = sum(CE * wts),
                    # CE_se
                    CE_se = sqrt(sum(CE_se ^ 2 * wts_sqr)),
                    CE_lo = NA_real_,
                    CE_hi = NA_real_, 
                    # uCE
                    uCE = sum(uCE * wts),
                    # uCE_se
                    uCE_se = sqrt(sum(uCE_se ^ 2 * wts_sqr)),
                    uCE_lo = NA_real_,
                    uCE_hi = NA_real_, 
                    # vCE
                    vCE = sum(vCE * wts),
                    # vCE_se
                    vCE_se = sqrt(sum(vCE_se ^ 2 * wts_sqr)),
                    vCE_lo = NA_real_, 
                    vCE_hi = NA_real_,
                    # wCE
                    wCE = sum(wCE * wts),
                    # wCE_se
                    wCE_se = sqrt(sum(wCE_se ^ 2 * wts_sqr)),
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
                    UCE = sum(UCE * wts)
                )
                # depositon:
                if (has_dep) {
                    out[, ':='(
                        # CE_Dep
                        CE_Dep = sum(CE_Dep * wts),
                        # CE_se_Dep
                        CE_se_Dep = sqrt(sum(CE_se_Dep ^ 2 * wts_sqr)),
                        # uCE_Dep
                        uCE_Dep = sum(uCE_Dep * wts),
                        # uCE_se_Dep
                        uCE_se_Dep = sqrt(sum(uCE_se_Dep ^ 2 * wts_sqr)),
                        # vCE_Dep
                        vCE_Dep = sum(vCE_Dep * wts),
                        # vCE_se_Dep
                        vCE_se_Dep = sqrt(sum(vCE_se_Dep ^ 2 * wts_sqr)),
                        # wCE_Dep
                        wCE_Dep = sum(wCE_Dep * wts),
                        # wCE_se_Dep
                        wCE_se_Dep = sqrt(sum(wCE_se_Dep ^ 2 * wts_sqr)),
                        # UCE_Dep
                        UCE_Dep = sum(UCE_Dep * wts)
                        )]
                }
                # all lo und hi nachrechnen... (TODO: aus by nehmen (N0m in out ausgeben)
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
            } else {
                cat('Fix me! .N != length(sensors_from)\n')
                browser()
            }
        }, by = .(rn, Source)]

        # get all point sensors
        Sensors <- attr(res, "ModelInput")$Sensors[attr(res, "ModelInput")$Sensors[, "Sensor Name"] %in% sensors_from, ]
        # old sensors
        Calc.Sensors <- procSensors(Sensors)$Calc.Sensors
        point_sensors <- Calc.Sensors[, "Point Sensor Name"]
        # new sensors
        NewSensors <- Sensors
        # IDs & Names
        NewSensors[, "Sensor ID"] <- paste(Sensors[, "Sensor Name"], Sensors[, "Sensor ID"], sep = "_")
        NewSensors[, "Sensor Name"] <- sensor_to
        NewCalc.Sensors <- procSensors(NewSensors)$Calc.Sensors
        new_point_sensors <- NewCalc.Sensors[, "Point Sensor Name"]
        names(new_point_sensors) <- point_sensors
        # fix CalcSteps
        CalcSteps_temp <- CalcSteps[Sensor %chin% sensors_from][, c('Sensor', 'Calc.Sensor') := {
            if (all(sensors_from %chin% Sensor)) {
                list(
                    Sensor = sensor_to,
                    Calc.Sensor = sapply(strsplit(Calc.Sensor, split = ','), function(x) {
                        paste(new_point_sensors[x], collapse = ',')
                    })
                )
            } else {
                NULL
            }
        }, by = .(rn, Source)]
        # Catalogs
        Catalogs_temp <- merge(Catalogs, CalcSteps_temp[, .(rn = unique(rn))], by = "rn")[PointSensor %chin% point_sensors]
        # replace Catalogs
        Catalogs_temp[, ':='(
            Sensor = sensor_to,
            PointSensor = new_point_sensors[PointSensor]
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
