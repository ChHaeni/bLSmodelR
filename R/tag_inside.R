# TODO:
#   - add argument to name output columns
#   - use better names for Sensors class (switchNames)
#       -> this must be switched right at beginning (in run_bls)
#       -> same for Sources -> name, x, y, id
#   - use tag_inside in calcCE.R ?
#   - use better classes than 'Sensors', 'Sources', 'TDcat' ('bls_sensors', etc.?)

#' tag touchdowns inside source areas
#'
#' tag touchdowns inside source areas. Adds columns
#'  "td_inside" and "source_names" to the TDcat object
#'
#' @param catalog A TDcat object, i.e., a bLSmodelR TD catalog of class TDcat.
#' @param sources An object of class "Sources".
#' @param origin Coordinates of the touchdown catalogs origin (sensor position) in 
#'  the coordinates of the "Sources" object.
#' @param tag_id Should polygon IDs be part of the source name? (Default: FALSE)
#' @return The provided TDcat object with two columns named "td_inside", indicating if
#'  touchdowns are inside the sources or not (TRUE/FALSE), and "source_names", indicating
#'  the name of the corresponding source where the touchdown hits inside.
tag_inside <- function(catalog, sources, origin = c(0, 0), tag_id = FALSE, 
    colname_inside = 'td_inside', colname_sources = 'source_names', name_outside = '') {

    if (!inherits(catalog, 'TDcat')) {
        stop('Argument "catalog" must be of class "TDcat"')
    }

    if (!inherits(sources, 'Sources')) {
        stop('Argument "sources" must be of class "Sources"')
    }
    
    if (!is.numeric(origin)) {
        if (!all(c('x-Coord (m)', 'y-Coord (m)') %in% names(origin))) {
            stop('Argument "origin" must be either a numeric vector of length 2 OR a ',
                ' matrix-like object containing a SINGLE ROW as well as COLUMN NAMES ',
                '"x-Coord (m)" and "y-Coord (m)"')
        }
        origin <- unlist(origin[, c('x-Coord (m)', 'y-Coord (m)'), drop = TRUE])
    }

    if (!is.numeric(origin) || length(origin) != 2) {
        stop('Argument "origin" must be either a numeric vector of length 2 OR a ',
            ' matrix-like object containing a SINGLE ROW as well as COLUMN NAMES ',
            '"x-Coord (m)" and "y-Coord (m)"')
    }

    if (!length(colname_inside) == 1 || !is.character(colname_inside)) {
        stop('Argument "colname_inside" must be a character vector of length == 1')
    }
    if (!length(colname_sources) == 1 || !is.character(colname_sources)) {
        stop('Argument "colname_sources" must be a character vector of length == 1')
    }

	catalog[, ':='(
        td_inside = FALSE,
        source_names = name_outside
        )]
	sources_relative <- data.table(sources)
	setnames(sources_relative, c("area", "x", "y", "pid"))
	sources_relative[, ":="(
        x = x - origin[[1]],
        y = y - origin[[2]]
        )]
	tag_bbox(catalog, sources_relative)
	catalog[, rn := .I]
	setkey(catalog, rn)

	if (catalog[, any(bbox_inside)]) {
        # assign outer bbox
        catalog[, bbox_outer := bbox_inside]
		# tag Inside Source
        if (tag_id) {
            tds_inside <- sources_relative[, {
                #reassing outer bbox and tag sub source
                catalog[, bbox_inside := bbox_outer]
                tag_bbox(catalog, .(x = x, y = y))
                cbind(
                    ID = catalog[(bbox_inside), rn], 
                    pnt.in.poly(catalog[(bbox_inside), cbind(x, y)], cbind(x, y))
                )
            }, by = .(area, pid)][pip == 1L, 
            paste0(paste(area, pid, sep = "/"), collapse = ",")
            , by = ID]	
        } else {
            tds_inside <- sources_relative[, {
                #reassing outer bbox and tag sub source
                catalog[, bbox_inside := bbox_outer]
                tag_bbox(catalog, .(x = x, y = y))
                cbind(
                    ID = catalog[(bbox_inside), rn], 
                    pnt.in.poly(catalog[(bbox_inside), cbind(x, y)], cbind(x, y))
                )
            }, by = .(area, pid)][pip == 1L, 
            paste0(area, collapse = ",")
            , by = ID]					
        }
        catalog[, bbox_outer := NULL]
        if (nrow(tds_inside) > 0) {
            catalog[tds_inside, ':='(
                td_inside = TRUE,
                source_names = V1
                )]
        }
	}
	catalog[, ":="(bbox_inside = NULL, rn = NULL)]

    setnames(catalog, c('td_inside', 'source_names'), c(colname_inside, colname_sources))

	invisible(catalog)
}


