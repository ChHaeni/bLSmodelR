
calc_fct <- function(res, ncores = NULL, SourceSplit = NULL, 
  dxy = 1, # in m
  tolerance = 1e-5, # digits for rounding
  checkArea = TRUE,
  memory_limit = NULL
  ){
  
  if(!is.null(ncores)) setDTthreads(ncores)

  # convert old versions 
  sres <- as.character(substitute(res))
  res <- copy(res)
  if(is.null(attr(res, "Version"))){
    warning(paste0("Object '", sres[min(length(sres), 2)], "' has not yet been converted to version 4.2+"))
    convert(res)
  }

  for(g in 1:20) gc()

  # get relevant sources
  uSous <- res[, unique(Source)]
  Sources <- attr(res, "ModelInput")$Sources
  SourcesList <- lapply(uSous, function(x){
    out <- Sources[Sources[,1] %in% x,]
    attr(out, "area") <- getArea(out)
    out
    })
  names(SourcesList) <- uSous

  # get relevant sensors
  uSens <- res[, unique(Sensor)]
  Sensors <- attr(res, "ModelInput")$Sensors
  C_sens_all <- procSensors(Sensors)$Calc.Sensors
  C_sens <- C_sens_all[C_sens_all[, "Sensor Name"] %in% uSens, ]

  # get relevant catalogs
  uRn <- res[, unique(rn)]
  All_cats <- Catalogs(res)[Sensor %in% uSens & rn %in% uRn]
  # prepare catalogs
  Cat_Path <- CatPath(res)
  All_cats[, CatNum := as.integer(as.factor(Cat.Name))]
  C_cats <- All_cats[,{
    .(
      CatNameList = list(setNames(Cat.Name, PointSensor)),
      CatNums = paste(sort(unique(CatNum)), collapse = ","),
      Subset_seed = seed[1]
      )
  }, by = .(rn, Sensor)]

  # split sources into sub areas
  if(is.null(SourceSplit)) {
    # split
    SourceSplit <- splitSource(Sources[Sources[, 1] %in% uSous, ],
      dxy = dxy, tolerance = tolerance, checkArea = checkArea)
  } else {
    for(i in seq_along(SourceSplit)){
      for(j in seq_along(SourceSplit[[i]])){
        # allocate memory
        SourceSplit[[i]][[j]]$polygons <- alloc.col(SourceSplit[[i]][[j]]$polygons)
      }
    }
  }
  
  # copy original data.table & add columns
  Calc <- merge(res[,.(
    rn,
    Sensor,
    Source,
    # Subset_seed,
    WD,
    SensorHeight,
    Ustar,
    L,
    Zo,
    sUu,
    sVu,
    bw,
    C0,
    kv,
    A,
    alpha,
    MaxFetch,
    N0,
    fct = NA_real_
    )], C_cats, by = c("rn", "Sensor"))

  if(!is.null(ncores) && ncores > 1){
    # showConnections(T)
    # getConnection(3)
    cl <- .makePSOCKcluster(ncores, memory_limit = memory_limit)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, {
      # library(data.table)
      # library(bLSmodelR)
      data.table::setDTthreads(1)
    })
    # call parallel
    InList <- lapply(seq_len(nrow(Calc)), function(x, y) y[x, ], y = as.data.frame(Calc))
    parallel::clusterExport(cl, c("SourceSplit", "SourcesList", "C_sens"), envir = environment())
    fct <- unlist(.clusterApplyLB(cl, InList, .fct, 
      C_Path = Cat_Path, parallel = TRUE))
    parallel::stopCluster(cl)
    on.exit()
  } else {

    for(g in 1:20) gc()

    # loop over rows
    fct <- .fct(Calc, SouSp = SourceSplit, SouLi = SourcesList, 
      c_sens = C_sens, C_Path = Cat_Path)
  }

  cbind(res, fct = fct)
}

# loop over rows
.fct <- function(DT, SouSp = SourceSplit, SouLi = SourcesList, c_sens = C_sens, 
  C_Path = Cat_Path, parallel = FALSE){
  if(parallel){
    setDT(DT)
    for(i in seq_along(SouSp)){
      for(j in seq_along(SouSp[[i]])){
        # allocate memory
        SouSp[[i]][[j]]$polygons <- alloc.col(SouSp[[i]][[j]]$polygons)
      }
    }
  }

  DT[,{
    if(!parallel){
      cat("*** row", .BY[[1]], "\n")
    }
    # get source
    iSource <- SouSp[[Source]]
    # add inside tag
    for(j in seq_along(iSource)){
      iSource[[j]]$polygons[, isInside := FALSE]
    }
    
    # list with point sensor names
    Point_Sensor <- split(names(CatNameList[[1]]), CatNameList[[1]])

    # loop over Catalogs
    for(Cat_Name in unique(CatNameList[[1]])){
      # get Catalog
      Cat0 <- readCatalog(file.path(C_Path, Cat_Name))
      # initialize
      initializeCatalog(.SD, Catalog = Cat0)
      # get cat N0
      Cat_N0 <- attr(Cat0, "N0")
      # take subset
      if (Cat_N0 > N0) {
        env <- globalenv()
        oseed <- env$.Random.seed
        set.seed(Subset_seed, kind = "L'Ecuyer-CMRG")      
        takeSub <- sample.int(Cat_N0, N0)
        if (is.null(oseed)) {
            rm(list = ".Random.seed", envir = env)
        } else {
            assign(".Random.seed", value = oseed, envir = env)
        }
        indexNew <- 1:N0
        names(indexNew) <- takeSub
        Cat0 <- Cat0[Traj_ID %in% takeSub, ]
        Cat0[, Traj_ID := indexNew[as.character(Traj_ID)]]
      }
      # rotate catalog
      rotateCatalog(Cat0, WD)
      
      # loop over point sensors
      for(ps in Point_Sensor[[Cat_Name]]){
        cat("|\t", ps, "\n")
        Sens_xy <- as.numeric(c_sens[c_sens[, "Point Sensor Name"] %in% ps, 
          c("x-Coord (m)", "y-Coord (m)")])
        # copy catalog
        Cat <- Cat0[, .(x, y, bbox_inside = TRUE, inside0 = FALSE)]
        # add Sensor x & y
        Cat[, ":="(
          x = x + Sens_xy[1],
          y = y + Sens_xy[2]
          )]
        # get IDs inside, initial inside has to be TRUE
        tag_bbox(Cat, apply(SouLi[[Source[1]]][,2:3], 2, range))
        Cat[, inside0 := bbox_inside]
        # initialize all_check
        all_check <- 0
        # loop over sources
        for(jSrc in seq_along(iSource)){
          # jSrc <- 1
          iSource[[jSrc]]$polygons[(!isInside), isInside := {
            # copy from original inside
            Cat[, bbox_inside := inside0]
            # tag near
            tag_bbox(Cat[(bbox_inside)], .(x, y))
            C_sub <- Cat[(bbox_inside), .(x, y)]
            any(as.logical(.Call("pip", C_sub[, x], C_sub[, y], C_sub[, .N], x, y, length(x), PACKAGE = "bLSmodelR")))
          }, by = .(tile, polygon)]
          all_check <- all_check + iSource[[jSrc]]$polygons[, as.integer(all(isInside))]
        }
        # stop if all tiles covered
        if(all_check == length(iSource)) break
      }
      # stop if all tiles covered
      if(all_check == length(iSource)) break
    }
    # calculate & return fct
    .(
      fct = sum(
        unlist(
          lapply(iSource, function(X) 
            X$polygons[(isInside), area[1], by = .(tile, polygon)][, sum(V1)]
          )
        )
      ) / attr(SouLi[[Source[1]]], "area")
    )
  }, by = .(.Row__ = seq_len(nrow(DT)))][order(.Row__), fct]
}


splitSource <- function(
  SourceIn,
  dxy = 1, # in m
  tolerance = 1e-5, # digits for rounding
  checkArea = TRUE
  ) {
  digits <- -floor(log10(tolerance))
  tol <- 10 ^ -digits

  setNumericRounding(2)

  # get source names
  uSourceNames <- unique(SourceIn[, 1])

  out <- lapply(uSourceNames, function(SourceName){
    # SourceName <- Source[1,1]
    Source <- SourceIn[SourceIn[, 1] %in% SourceName, ]
    # start lapply over polyID
    lapply(unique(Source[,4]), function(polyID){
      SubSource <- Source[Source[,4] %in% polyID,]
      xyMeans <- colMeans(SubSource[,2:3])
      SubSource[,2:3] <- round(sweep(SubSource[, 2:3], 2, xyMeans), digits - 1)
      Poly <- setNames(SubSource[,2:3], c("x", "y"))
      Poly_extr <- attr(procSources(SubSource), "SourceList")[[SourceName]][, 1:2]
      Poly_lines <- lapply(seq_len(nrow(Poly_extr) - 1), function(x) {
        list(
          p1 = as.numeric(Poly_extr[x, ]),
          p2 = as.numeric(Poly_extr[x + 1, ])
          )
      })
      # # possibly problems with values close to 0?
      # Poly_lines_dir <- lapply(Poly_lines, function(x){
      #   list(x = sign(x$p2[1] - x$p1[1]), y = sign(x$p2[2] - x$p1[2]))
      # })

      ##### Ablauf:
      # browser()
      # get BB (extended)
      BBoxOrig <- list(
        x = trunc(range(Poly$x)) + c(-1, 1), 
        y = trunc(range(Poly$y)) + c(-1, 1)
        # x = range(Poly$x), 
        # y = range(Poly$y)
        )
      BBdim <- lapply(BBoxOrig, diff)
      BBdim_ext <- lapply(BBdim, function(x) ceiling(x / dxy) * dxy)
      BBadd <- mapply("-", BBdim_ext, BBdim) / 2
      BBox <- mapply(function(x, y) x + c(-1, 1) * y, x = BBoxOrig, y = BBadd, SIMPLIFY = FALSE)

      # liste mit matrix mit dxy -> grid
      # x & y vectors
      x_seq <- seq(BBox$x[1], BBox$x[2], by = dxy)
      y_seq <- seq(BBox$y[1], BBox$y[2], by = dxy)

      # create matrix M
      dimx <- length(x_seq) - 1
      dimy <- length(y_seq) - 1
      M <- matrix(seq_len(dimx * dimy), ncol = dimy, nrow = dimx)

      ####### inside
      # create grids Gxy & G
      Gxy <- as.matrix(expand.grid(row = seq_along(x_seq), col = seq_along(y_seq), KEEP.OUT.ATTRS = FALSE))
      G <- as.matrix(expand.grid(x = x_seq, y = y_seq, KEEP.OUT.ATTRS = FALSE))

      # - check Eckpunkte inside
      corner_inside <- as.logical(.Call("pip", G[, 1], G[, 2], nrow(G), Poly[, 1], 
        Poly[, 2], nrow(Poly), PACKAGE = "bLSmodelR"))
      corner_inside_index <- which(corner_inside)
      Gxy_inside <- Gxy[corner_inside,, drop = FALSE]

      # get tiles in M
      M_xy <- apply(Gxy_inside, 1, function(x){
        c(x[1] + c(-1, 0), x[2] + c(-1, 0))
      })
      M_xy[M_xy == 0] <- NA
      M_tiles <- apply(M_xy, 2, function(x){
        c(
          (x[3] - 1) * dimx + x[1],
          (x[3] - 1) * dimx + x[2],
          (x[4] - 1) * dimx + x[1],
          (x[4] - 1) * dimx + x[2]
          )
      })
      # count corners inside per tile
      Tiles_count <- table(M_tiles)

     ####### intersecting
      # get all points with x_i | y_j intersect mit source_poly
      vert <- lapply(x_seq, function(px) list(p1 = c(px, BBox$y[1]), p2 = c(px, BBox$y[2])))
      horiz <- lapply(y_seq, function(py) list(p1 = c(BBox$x[1], py), p2 = c(BBox$x[2], py)))
      
      # check vertical grid
      Vcheck <- do.call("rbind", 
        lapply(seq_along(vert), function(v) do.call("rbind", 
          lapply(seq_along(Poly_lines), function(i){
            cbind(inters(Poly_lines[[i]]$p1, Poly_lines[[i]]$p2, vert[[v]]$p1, vert[[v]]$p2, tol), 
              Poly_lines_entry = i, Line_row = v, Line_col = NA)
      }))))
      # check horizontal grid
      Hcheck <- do.call("rbind", 
        lapply(seq_along(horiz), function(v) do.call("rbind",
          lapply(seq_along(Poly_lines), function(i){
            cbind(inters(Poly_lines[[i]]$p1, Poly_lines[[i]]$p2, horiz[[v]]$p1, horiz[[v]]$p2, tol), 
              Poly_lines_entry = i, Line_row = NA, Line_col = v)
      }))))

      # bind together
      All_points <- rbind(
        Vcheck[Vcheck[,3] == 1, ], 
        Hcheck[Hcheck[,3] == 1, ]
        )

      # inters(Poly_lines[[1]]$p1, Poly_lines[[1]]$p2, horiz[[18]]$p1, horiz[[18]]$p2, tol)

      # find segment lines by x & y -> row, col
      Line_rows <- All_points[,"Line_row"]
      isna_row <- is.na(Line_rows)
      Line_rows[isna_row] <- findInterval(All_points[isna_row, 1], x_seq, rightmost.closed = TRUE)
      Line_cols <- All_points[,"Line_col"]
      isna_col <- is.na(Line_cols)
      Line_cols[isna_col] <- findInterval(All_points[isna_col, 2], y_seq, rightmost.closed = TRUE)

      # check exact matching
      Outer_rows <- abs(outer(All_points[, 1], x_seq, "-")) < tol
      Exact_rows <- as.logical(rowSums(Outer_rows) > 0)
      Outer_cols <- abs(outer(All_points[, 2], y_seq, "-")) < tol
      Exact_cols <- as.logical(rowSums(Outer_cols) > 0)
      # check for exact matches in foundIntervals
      if(any(Exact_rows[isna_row])){
        Line_rows[Exact_rows & isna_row] <- apply(Outer_rows[Exact_rows & isna_row, , drop = FALSE], 1, which)
      }
      if(any(Exact_cols[isna_col])){
        Line_cols[Exact_cols & isna_col] <- apply(Outer_cols[Exact_cols & isna_col, , drop = FALSE], 1, which)
      }

      # check exact corner
      Exact_corners <- Exact_rows & Exact_cols
      
      # round with digits
      All_points[, "x"] <- round(All_points[, "x"], digits - 1)
      All_points[, "y"] <- round(All_points[, "y"], digits - 1)

      # find Tiles to check: lines -> 2 Tiles
      Border_x <- rbind(
        data.table(
          row = Line_rows[Exact_rows & !Exact_corners],
          col = Line_cols[Exact_rows & !Exact_corners],
          Px = All_points[Exact_rows & !Exact_corners, 1],
          Py = All_points[Exact_rows & !Exact_corners, 2],
          Poly_lines_entry = All_points[Exact_rows & !Exact_corners, 4]
          ),
        data.table(
          row = Line_rows[Exact_rows & !Exact_corners] - 1,
          col = Line_cols[Exact_rows & !Exact_corners],
          Px = All_points[Exact_rows & !Exact_corners, 1],
          Py = All_points[Exact_rows & !Exact_corners, 2],
          Poly_lines_entry = All_points[Exact_rows & !Exact_corners, 4]
          ))
      Border_y <- rbind(
        data.table(
          row = Line_rows[Exact_cols & !Exact_corners],
          col = Line_cols[Exact_cols & !Exact_corners],
          Px = All_points[Exact_cols & !Exact_corners, 1],
          Py = All_points[Exact_cols & !Exact_corners, 2],
          Poly_lines_entry = All_points[Exact_cols & !Exact_corners, 4]
          ),
        data.table(
          row = Line_rows[Exact_cols & !Exact_corners],
          col = Line_cols[Exact_cols & !Exact_corners] - 1,
          Px = All_points[Exact_cols & !Exact_corners, 1],
          Py = All_points[Exact_cols & !Exact_corners, 2],
          Poly_lines_entry = All_points[Exact_cols & !Exact_corners, 4]
          ))

      # find Tiles to check: corners -> 4 Tiles
      Border_xy <- rbind(
        data.table(
          row = Line_rows[Exact_corners],
          col = Line_cols[Exact_corners],
          Px = All_points[Exact_corners, 1],
          Py = All_points[Exact_corners, 2],
          Poly_lines_entry = All_points[Exact_corners, 4]
          ),
        data.table(
          row = Line_rows[Exact_corners] - 1,
          col = Line_cols[Exact_corners],
          Px = All_points[Exact_corners, 1],
          Py = All_points[Exact_corners, 2],
          Poly_lines_entry = All_points[Exact_corners, 4]
          ),
        data.table(
          row = Line_rows[Exact_corners],
          col = Line_cols[Exact_corners] - 1,
          Px = All_points[Exact_corners, 1],
          Py = All_points[Exact_corners, 2],
          Poly_lines_entry = All_points[Exact_corners, 4]
          ),
        data.table(
          row = Line_rows[Exact_corners] - 1,
          col = Line_cols[Exact_corners] - 1,
          Px = All_points[Exact_corners, 1],
          Py = All_points[Exact_corners, 2],
          Poly_lines_entry = All_points[Exact_corners, 4]
          )
        )

      # gather all tiles
      Border_all <- rbind(Border_x,Border_y,Border_xy)

      # remove borders
      Border <- Border_all[
        row > 0 &
        col > 0 &
        row <= dimx &
        col <= dimy,]
      # add Tile number
      Border[, tile := (Border$col - 1) * dimx + Border$row]

      # get number of corners inside source
      Border[, nc := 0L]
      Border[as.character(tile) %chin% names(Tiles_count), nc := {
            t_c <- as.character(tile)
            Tiles_count[t_c[t_c %in% names(Tiles_count)]]
            }]


      # # -> was falls mehrere polygons? (eher unwahrsch.)
      # # einbauen: if(count(Poly_lines_entry) > 2) -> mehrere Polygons
      # if(polyID == 2)browser()

      # check & get polygons
      All <- Border[, {

        # borders
        x_le <- x_seq[row[1]]
        x_ri <- x_seq[row[1] + 1]
        y_lo <- y_seq[col[1]]
        y_hi <- y_seq[col[1] + 1]

        # any points inside?
        pts_inside <- suppressWarnings(do.call(rbind, lapply(seq_len(nrow(Poly)), function(i){
            if(Poly[i, 1] >= x_le & Poly[i, 1] <= x_ri & 
              Poly[i, 2] >= y_lo & Poly[i, 2] <= y_hi ){
              data.table(Poly[i, ], label = paste0("L",c((i - 2) %% nrow(Poly) + 1, i)))
            } else {
              NULL
            }
          })))

        # get corners inside
        Corner <- as.data.table(G[corner_inside_index[which(colSums(M_tiles == .BY[[1]]) > 0)], , drop = FALSE])

        # get labels of corners
        Corner[, label := {
          LR <- c("L", "R")[(abs(x - x_ri) < tol) + 1]
          BT <- c("B", "T")[(abs(y - y_hi) < tol) + 1]
          paste0("C", sapply(paste0(LR,BT), switch,
            "LB" = 1,
            "RB" = 2,
            "RT" = 3,
            "LT" = 4
            ))
        }]

        # browser()

        # cat(.BY$tile,"\n")
        # if(.BY$tile == 14) browser() # -> Corner of Source Area!
        # if(.BY$tile == 18) browser() # -> sprintf of 0 results in -0.00 ???!!!

        # check if line corner is exactly on corner (klappt das so, falls beide lines entlang dem edge gehen???)
        fmt <- paste0("%1.", digits,"f")
        fmt_xy <- paste(fmt, fmt, sep = "::")
        crn_xy <- Corner[, sprintf(fmt_xy, abs(x) * sign(x), abs(y) * sign(y))]
        sd_xy <- sprintf(fmt_xy, abs(Px) * sign(Px), abs(Py) * sign(Py))
        if(any(sd_dup <- duplicated(sd_xy))){
          sd_check <- sd_xy[sd_dup]
          sd_remove <- sd_check[sd_check %in% crn_xy]
          if(length(sd_remove) > 0){
            # remove lines on corner
            SD0 <- .SD[!(sd_xy %in% sd_remove), 
                .(x = Px, y = Py, label = paste0("L", Poly_lines_entry))]
            sd_xy <- sd_xy[!(sd_xy %in% sd_remove)]
          } else {
            # remove corners on line
            SD0 <- .SD[, .(x = Px, y = Py, label = paste0("L", Poly_lines_entry))]
          }
        } else if(.N > 1) {
          # remove corners on line
          SD0 <- .SD[, .(x = Px, y = Py, label = paste0("L", Poly_lines_entry))]
        } else {
          # corner lies ecaxtly on line corner
          SD0 <- NULL
          sd_xy <- NULL
        }

        # check if corner falls on line or on pnts inside
        if(is.null(pts_inside)){
          if(any(crn_xy %in% sd_xy)) {
            Corner0 <- Corner[!(crn_xy %in% sd_xy)]
          } else {
            Corner0 <- Corner
          }
        } else {
          inside_xy <- pts_inside[, sprintf(fmt_xy, abs(x) * sign(x), abs(y) * sign(y))]
          if(any(crn_xy %in% sd_xy, crn_xy %in% inside_xy)) {
            Corner0 <- Corner[!(crn_xy %in% sd_xy) & !(crn_xy %in% inside_xy)]
          } else {
            Corner0 <- Corner
          }      
        }

        # rbind (unique nicht mehr noetig?)
        DT <- unique(rbind(
          Corner0,
          # Corner,
          # if(!is.null(SD0)) SD0[.N:1, ],
          SD0,
          pts_inside
          ))

        # plot(c(x_le, x_ri), c(y_lo, y_hi), type = "n")
        # rect(x_le, y_lo, x_ri, y_hi, lwd = 2, border = "darkgrey")
        # if(!is.null(pts_inside)) pts_inside[, points(x, y, col = "darkgreen")]
        # Corner[, points(x, y, pch = 20, cex = 1.5)]
        # PLE <- unique(Poly_lines_entry)
        # for(ple in seq_along(PLE)){
        #   points(Px[Poly_lines_entry == PLE[ple]], Py[Poly_lines_entry == PLE[ple]], type = "b", col = "red")
        # }

        # if(.BY$tile == 9721) browser()
        # if(polyID == 2 && .BY$tile == 187) browser()

        # go around corners
        if(DT[, .N > 0]){
          DT[, ":="(
            included = FALSE,
            rn = 1:.N,
            polygon = NA_character_,
            polygon_order = NA_integer_,
            edge = checkEdge(x, y, y_lo, x_ri, y_hi, x_le, tol),
            corners = nc[1]
            ) ]
          poly_order <- poly_number <- 1L
          # get first point
          pnt <- DT[1,]
          # set polygon values in DT
          DT[1, ":="(
            included = TRUE,
            polygon = letters[poly_number],
            polygon_order = poly_order
            )]
          while(DT[,any(!included)]){
            # get point label
            pnt_lbl <- pnt[, label]
            # check if corner or line
            if(grepl("C", pnt_lbl)){
              # corner
              if(pnt_lbl %in% c("C4", "C2")){
                # check matching x & best y
                pnt <- DT[!included & abs(x - pnt[, x]) < tol][
                  which.min(abs(y - pnt[, y]))]
              } else {
                # check matching y & best x
                pnt <- DT[!included & abs(y - pnt[, y]) < tol][
                  which.min(abs(x - pnt[, x]))]
              }
            } else {
              # line
              if(DT[!(included), pnt_lbl %chin% label]){
                # take next point with label
                pnt <- DT[!included & label %in% pnt_lbl]
              } else {
                # search for next point on corresponding edge
                if(pnt[, edge == 5]){
                  # point inside -> check for point repetition
                  pnt <- DT[!included & abs(x - pnt[,x]) < tol & abs(y - pnt[,y]) < tol]
                } else {
                  # check side for next point
                  Side <- pnt[, edge]
                  if(Side %in% c(2, 4)){
                    # check matching x & best y
                    pnt <- DT[!included & abs(x - pnt[, x]) < tol][
                      which.min(abs(y - pnt[, y]))]
                  } else {
                    # check matching y & best x
                    pnt <- DT[!included & abs(y - pnt[, y]) < tol][
                      which.min(abs(x - pnt[, x]))]
                  }
                }
              }
            }
            if(nrow(pnt)){
              # update poly_order
              poly_order <- poly_order + 1
              # set polygon values in DT
              DT[rn == pnt[, rn], ":="(
                included = TRUE,
                polygon = letters[poly_number],
                polygon_order = poly_order
                )]
            } else {
              # update poly_* for next polygon
              poly_order <- 1
              poly_number <- poly_number + 1
              pnt <- DT[!(included)][1,]     
              DT[rn == pnt[, rn], ":="(
                included = TRUE,
                polygon = letters[poly_number],
                polygon_order = poly_order
                )]
            }
          }


          # DT[order(polygon_order)][,{
          #   polygon(x, y, 
          #     col = c("#D63636", "#2E64EE", "green", "orange", "purple")[.BY$corners + 1])
          # }, by = .(corners, polygon)]

          # DT[order(polygon_order), plot(x, y, type = "b")]
          DT[order(polygon_order), ][, area := {
            ind <- c(.N, seq(.N - 1))
            abs(sum(x * y[ind] - x[ind] * y)) / 2          
            }, by = polygon]
        } else {
          NULL
        }
      }, by = tile]

      # correct for xyMeans
      x_seq <- x_seq + xyMeans[1]
      y_seq <- y_seq + xyMeans[2]
      Poly$x <- Poly$x + xyMeans[1]
      Poly$y <- Poly$y + xyMeans[2]

      # output
      out <- list(
        x_seq = x_seq,
        y_seq = y_seq,
        polygons = rbind(
          All[, .(tile, polygon,
            x = x + xyMeans[1], 
            y = y + xyMeans[2], 
            corners, area)],
          data.table(tile = as.integer(names(Tiles_count[Tiles_count == 4])) %w/o% All[,unique(tile)])[, {
          col <- ceiling(tile / dimx)
          row <- tile - (col - 1) * dimx
          .(
            polygon = "a",
            x = x_seq[row + c(0, 1, 1, 0)],
            y = y_seq[col + c(0, 0, 1, 1)],
            corners = -1L,
            area = dxy * dxy
          )}, by = tile]),
        dxy = dxy,
        tolerance = tolerance,
        digits = digits,
        tol = tol,
        SourcePoly = Poly
        )

      if(checkArea && abs(getArea(SubSource) - out$polygons[,area[1], by = .(tile, polygon)][,sum(V1)]) > sqrt(tol)){
        x11()
        # plotM(out)
        plot(Poly$x[c(1:nrow(Poly),1)], Poly$y[c(1:nrow(Poly),1)], type = "l", lwd = 2)
          # ,xlim = c(-1, 1) + 711620, ylim = c(-1, 1) + 260850)
        plotPolys(out)
        stop("polygons area and total area differ by more than sqrt(tol)! try to set a better dxy..")
      } else {
        out
      }

    })
    # end lapply over polygon ID
  })
  # end lapply over sources
  setNames(out, uSourceNames)
}

inters <- function(l1_p1, l1_p2, l2_p1, l2_p2, tol = 1e-5){
  # rewrite as vectors
  # -> https://demonstrations.wolfram.com/IntersectionOfTwoLinesUsingVectors/
  # browser()
  VZ <- l1_p2 - l1_p1
  UW <- l2_p2 - l2_p1
  mvz <- sqrt(sum(VZ ^ 2))
  muw <- sqrt(sum(UW ^ 2))
  beta <- acos(sum(VZ * UW) / (mvz * muw))
  UV <- l2_p1 - l1_p1
  muv <- sqrt(sum(UV ^ 2))
  alpha <- acos(sum(UV * UW) / (muv * muw))

  out <- l1_p1 + muv * sin(alpha) / sin(beta) * VZ / mvz
  
  return(cbind(
    x = out[1],
    y = out[2],
    inside =
      (
        (
          (out[1] - l1_p1[1]) >= -tol &&
          (out[1] - l1_p2[1]) <= tol
        ) || (
          (out[1] - l1_p1[1]) <= tol &&
          (out[1] - l1_p2[1]) >= -tol
        )
      ) && (
        (
          (out[1] - l2_p1[1]) >= -tol &&
          (out[1] - l2_p2[1]) <= tol
        ) || (
          (out[1] - l2_p1[1]) <= tol &&
          (out[1] - l2_p2[1]) >= -tol
        )
      ) && (
        (
          (out[2] - l1_p1[2]) >= -tol &&
          (out[2] - l1_p2[2]) <= tol
        ) || (
          (out[2] - l1_p1[2]) <= tol &&
          (out[2] - l1_p2[2]) >= -tol
        )
      ) && (
        (
          (out[2] - l2_p1[2]) >= -tol &&
          (out[2] - l2_p2[2]) <= tol
        ) || (
          (out[2] - l2_p1[2]) <= tol &&
          (out[2] - l2_p2[2]) >= -tol
        )
      ) 
    ))
}

checkEdge <- function(x, y, y_lo, x_ri, y_hi, x_le, tol = 1e-5){
  out <- rep(1, length(x))
  out[abs(x - x_le) < tol] <- out[abs(x - x_le) < tol] * 11
  out[abs(y - y_hi) < tol] <- out[abs(y - y_hi) < tol] * 7
  out[abs(x - x_ri) < tol] <- out[abs(x - x_ri) < tol] * 5
  out[abs(y - y_lo) < tol] <- out[abs(y - y_lo) < tol] * 3
  sapply(as.character(out), switch,
    "3" = 1, "33" = 1,
    "5" = 2, "15" = 2,
    "7" = 3, "35" = 3,
    "11" = 4, "77" = 4, 5)
}

# plotM <- function(X, xlim = NULL, ylim = NULL, cex.text = 0.35){
#   plot(range(X$x_seq), range(X$y_seq),type = "n", xlim = xlim, ylim = ylim)
#   lines(X$SourcePoly$x[c(1:nrow(X$SourcePoly),1)],X$SourcePoly$y[c(1:nrow(X$SourcePoly),1)], lwd = 2)
#   abline(h = X$y_seq)
#   abline(h = X$y_seq[c(1, length(X$y_seq))], lwd = 2)
#   abline(v = X$x_seq)
#   abline(v = X$x_seq[c(1, length(X$x_seq))], lwd = 2)
#   X$polygons[,{
#     polygon(x, y, 
#       col = c("lightgrey", "#D63636", "#2E64EE", "green", "orange", "purple")[.BY$corners + 2])
#   }, by = .(corners, tile, polygon)]
#   Tgrid <- expand.grid(X$x_seq[-1] - X$dxy/2, X$y_seq[-1] - X$dxy/2)
#   text(Tgrid[,1], Tgrid[,2], 1:nrow(Tgrid), cex = cex.text)
# }

plotPolys <- function(X, xlim = NULL, ylim = NULL, cex.text = 0.35,
  col = c("lightgrey", "#D63636", "#2E64EE", "green", "orange", "purple"),
  tileNumbers = TRUE){
  if(length(col) == 1) col <- rep(col , 6)
  X$polygons[,{
    polygon(x, y, 
      col = col[.BY$corners + 2])
  }, by = .(corners, tile, polygon)]
  if(tileNumbers){
    Tgrid <- expand.grid(X$x_seq[-1] - X$dxy/2, X$y_seq[-1] - X$dxy/2)
    text(Tgrid[,1], Tgrid[,2], 1:nrow(Tgrid), cex = cex.text)
  }
}
