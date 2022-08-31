combineSources <- function(res, comb_list = NULL, weight_units = c("m/A/t", "m/t"), add = FALSE){

  # convert old versions 
  sres <- as.character(substitute(res))
  res <- copy(res)
  setDT(res)
  if(is.null(attr(res, "Version"))){
    warning(paste0("Object '", sres[min(length(sres), 2)], "' has not yet been converted to version 4.2+"))
    convert(res)
  }


  if(is.null(comb_list)){
    uSources <- res[, unique(Source)]
    comb_list <- setNames(
      list(
        setNames(rep(1, length(uSources)), uSources)
        ), "all_sources_combined")
  } else {
    # check weights
    comb_list <- lapply(comb_list, function(x){
      if(is.character(x)){
        setNames(rep(1, length(x)), x)
      } else {
        if(weight_units[1] %in% "m/t"){
          source_areas <- unique(res, by = "Source")[Source %in% names(x), setNames(SourceArea, Source)]
          x / source_areas[names(x)]
        } else {
          x
        }
      }
    })
  }

  # check duplicates
  if(any(duplicated(c(names(comb_list), res[, unique(Source)])))) stop("please provide new unique names for combined sources")

  # do fct values exist?
  fct_exists <- "fct" %in% names(res)

  outlist <- setNames(vector("list", length(comb_list)), names(comb_list))
  for(source_to in names(comb_list)){
    # source_to <- names(comb_list)[1]
    sources <- names(comb_list[[source_to]])
    cat('Combining sources', paste(sources, collapse = ' + '), 'to', source_to, '\n')
    prior_wts <- comb_list[[source_to]]
    outlist[[source_to]] <- res[Source %in% sources, {
      if(.N == length(sources)){
        CE_wts <- avgSourcesWeights(prior_wts[Source], SourceArea)
        out <- data.table(
          # Source
          Source = source_to,
          # SensorHeight
          SensorHeight = SensorHeight[1],
          # SourceArea
          SourceArea = sum(SourceArea),
          # CE
          CE = avgCE_sources(CE, CE_wts),
          # CE_se
          CE_se = avgCE_sources_se(CE_se, CE_wts),
          CE_lo = NA_real_,
          CE_hi = NA_real_, 
          # uCE
          uCE = avgCE_sources(uCE, CE_wts),
          # uCE_se
          uCE_se = avgCE_sources_se(uCE_se, CE_wts),
          uCE_lo = NA_real_,
          uCE_hi = NA_real_, 
          # vCE
          vCE = avgCE_sources(vCE, CE_wts),
          # vCE_se
          vCE_se = avgCE_sources_se(vCE_se, CE_wts),
          vCE_lo = NA_real_, 
          vCE_hi = NA_real_,
          # wCE
          wCE = avgCE_sources(wCE, CE_wts),
          # wCE_se
          wCE_se = avgCE_sources_se(wCE_se, CE_wts),
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
          UCE = avgCE_sources(UCE, CE_wts),
          # fct
          fct = if(fct_exists) sum(fct) else NULL
          )
        # depositon:
        if (inherits(result, 'deposition')) {
            out[, ':='(
                  # CE_Dep
                  CE_Dep = avgCE_sources(CE_Dep, CE_wts),
                  # CE_se_Dep
                  CE_se_Dep = avgCE_sources_se(CE_se_Dep, CE_wts),
                  # uCE_Dep
                  uCE_Dep = avgCE_sources(uCE_Dep, CE_wts),
                  # uCE_se_Dep
                  uCE_se_Dep = avgCE_sources_se(uCE_se_Dep, CE_wts),
                  # vCE_Dep
                  vCE_Dep = avgCE_sources(vCE_Dep, CE_wts),
                  # vCE_se_Dep
                  vCE_se_Dep = avgCE_sources_se(vCE_se_Dep, CE_wts),
                  # wCE_Dep
                  wCE_Dep = avgCE_sources(wCE_Dep, CE_wts),
                  # wCE_se_Dep
                  wCE_se_Dep = avgCE_sources_se(wCE_se_Dep, CE_wts),
                  # UCE_Dep
                  UCE_Dep = avgCE_sources(UCE_Dep, CE_wts),
                )]
        }
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
        # cbind all missing columns!
        cbind(
          out,
          .SD[1, names(.SD) %w/o% names(out), with = FALSE]
          )
      }
    }, by = .(rn, Sensor)]
  }

  # adjust attributes part 1
  CalcSteps <- copy(attr(res, "CalcSteps"))
  ModelInput <- attr(res, "ModelInput")

  # prepare out
  out <- rbindlist(outlist)
  setcolorder(out, names(res))

  if(add){
    # add sources to Source
    for(source_to in names(comb_list)){
      # source_to <- names(comb_list)[1]
      pattern <- paste0("\\b", names(comb_list[[source_to]]), "\\b")
      # get rows
      ind <- CalcSteps[, do.call(mapply, c(FUN = all, lapply(pattern, function(x) grepl(x, Source))))]
      # CalcSteps
      CalcSteps[ind, Source := paste(Source, source_to, sep = ",")]
      # Sources
      Srcs <- ModelInput$Sources[ModelInput$Sources[, 1] %in% names(comb_list[[source_to]]),]
      Srcs[, 4] <- as.integer(as.factor(Srcs[, 1]))
      Srcs[, 1] <- source_to
      ModelInput$Sources <- rbind(
        ModelInput$Sources,
        Srcs
        )
    }
    out <- rbind(
      res,
      out
      )
  } else {
    CalcSteps[, Source__new := ""]
    # add sources to Source
    ModelInput_Sources <- NULL
    for(source_to in names(comb_list)){
      # source_to <- names(comb_list)[1]
      pattern <- paste0("\\b", names(comb_list[[source_to]]), "\\b")
      # get rows
      ind <- CalcSteps[, do.call(mapply, c(FUN = all, lapply(pattern, function(x) grepl(x, Source))))]
      # CalcSteps
      CalcSteps[ind, Source__new := paste(Source__new, source_to, sep = ",")]
      # Sources
      Srcs <- ModelInput$Sources[ModelInput$Sources[, 1] %in% names(comb_list[[source_to]]),]
      Srcs[, 4] <- as.integer(as.factor(Srcs[, 1]))
      Srcs[, 1] <- source_to
      ModelInput_Sources <- rbind(ModelInput_Sources, Srcs)
    }
    ModelInput$Sources <- ModelInput_Sources
    # get relevant rows and rewrite Source
    CalcSteps <- CalcSteps[Source__new != ""][, ":="(
      Source = Source__new,
      Source__new = NULL
      )]
    # remove first comma
    CalcSteps[, Source := sub("^[,]", "", Source)]
  }

  # adjust attributes part 2
  for(att in names(attributes(res)) %w/o% c("row.names", 
    ".internal.selfref", "names", "index", "CalcSteps", "ModelInput")){
    setattr(out, att, attr(res, att))
  }

  setkey(CalcSteps, rn, Sensor)

  setattr(out, "CalcSteps", CalcSteps)
  setattr(out, "ModelInput", ModelInput)
  setattr(out, "combined_sources_list", comb_list)
  
  setorder(out, rn, Sensor, Source)
  setindex(out, Source)

  out
}


avgSourcesWeights <- function(prior_weights, source_areas){
  R <- outer(prior_weights, prior_weights, "/")
  sum(source_areas) / as.numeric(t(R) %*% source_areas)
}

avgCE_sources <- function(ce, weights){
  sum(weights * ce)
}
avgCE_sources_se <- function(ce_se, weights){
  sqrt(sum(weights ^ 2 * ce_se ^ 2))
}

