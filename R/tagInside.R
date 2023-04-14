# TODO:
#   - change to consistent tagging (character instead of character/logical)
#       -> careful rattail!
#   - rename function 
#       -> rename all existing tagInside() calls
#   - better arguments/argument names
#   - use better names for Sensors class (switchNames)
#       -> this must be switched right at beginning (in run_bls)
#       -> same for Sources -> name, x, y, id
#   - use tag_inside in calcCE.R ?
#   - use better classes than 'Sensors', 'Sources', 'TDcat' ('bls_sensors', etc.?)

#' tag touchdowns inside source areas
#'
#' tag touchdowns inside source areas. Adds columns
#'  "tag_inside" and "source_names" to the TDcat object
#'
#' @param catalog A TDcat object, i.e., a bLSmodelR TD catalog of class TDcat.
#' @param sources An object of class "Sources".
#' @param origin Coordinates of the touchdown catalogs origin (sensor position) in 
#'  the coordinates of the "Sources" object.
#' @param tag_id Should polygon IDs be part of the source name? (Default: FALSE)
#' @return The provided TDcat object with two columns named "tag_inside", indicating if
#'  touchdowns are inside the sources or not (TRUE/FALSE), and "source_names", indicating
#'  the name of the corresponding source where the touchdown hits inside.
tagInside <- function(catalog, sources, origin = c(0, 0),
    tag_id = FALSE) {

    if (!inherits(catalog, 'TDcat')) {
        stop('Argument "catalog" must be of class "TDcat"')
    }

    if (!inherits(sources, 'Sources')) {
        stop('Argument "sources" must be of class "Sources"')
    }
    
    if (inherits(origin, 'Sensors')) {
        origin <- origin[, c('x-Coord (m)', 'y-Coord (m)'), drop = TRUE]
    }

    if (!is.numeric(origin) || length(origin) != 2) {
        stop('Argument "origin" must be either a numeric vector of length 2 or a single',
            ' point sensor (class "Sensors")')
    }

	catalog[, ':='(
        tag_inside = FALSE,
        source_names = ''
        )]
	sources_relative <- data.table(sources)
	setnames(sources_relative, c("area", "x", "y", "pid"))
	sources_relative[, ":="(
        x = x - origin[1],
        y = y - origin[2]
        )]
	tag_bbox(catalog, sources_relative)
	catalog[, rn := .I]
	setkey(catalog, rn)

	if (catalog[, any(inside)]) {
		# tag Inside Source
        if (tag_id) {
            td_inside <- sources_relative[, {
                cbind(
                    ID = catalog[(inside), rn], 
                    pnt.in.poly(catalog[(inside), cbind(x, y)], cbind(x, y))
                )
            }, by = .(area, pid)][pip == 1L, 
            paste0(paste(area, pid, sep = "/"), collapse = ",")
            , by = ID]	
        } else {
            td_inside <- sources_relative[, {
                cbind(
                    ID = catalog[(inside), rn], 
                    pnt.in.poly(catalog[(inside), cbind(x, y)], cbind(x, y))
                )
            }, by = .(area, pid)][pip == 1L, 
            paste0(area, collapse = ",")
            , by = ID]					
        }
        catalog[td_inside, ':='(
            tag_inside = TRUE,
            source_names = V1
            )]
	}
	catalog[, ":="(inside = NULL, rn = NULL)]

	invisible(catalog)
}


