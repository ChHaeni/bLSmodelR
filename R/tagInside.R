# TODO:
#   - change to consistent tagging (character instead of character/logical)
#       -> careful rattail!
#   - rename function 
#       -> rename all existing tagInside() calls
#   - better arguments/argument names
#   - use better names for Sensors class (simpleNames)
#       -> this must be switched right at beginning (in run_bls)
#   - use tag_inside in calcCE.R ?
#   - use better classes than 'Sensors', 'Sources', 'TDcat' ('bls_sensors', etc.?)
tagInside <- function(catalog, tag_source, origin = c(0, 0),
    tagBySourceName = tagPoly, tagPoly = FALSE) {

    if (!inherits(catalog, 'TDcat')) {
        stop('Argument "catalog" must be of class "TDcat"')
    }

    if (!inherits(tag_source, 'Sources')) {
        stop('Argument "tag_source" must be of class "Sources"')
    }
    
    if (inherits(origin, 'Sensors')) {
        origin <- origin[, c('x-Coord (m)', 'y-Coord (m)'), drop = TRUE]
    }

    if (!is.numeric(origin) || length(origin) != 2) {
        stop('Argument "origin" must be either a numeric vector of length 2 or a single',
            ' point sensor (class "Sensors")')
    }

	catalog[, tagInside := 
        if (tagBySourceName) {
            rep("", nrow(catalog))
        } else {
            rep(FALSE, nrow(catalog))
        }
    ]
	source_relative <- data.table(tag_source)
	setnames(source_relative, c("area", "x", "y", "pid"))
	source_relative[, ":="(
        x = x - origin[1],
        y = y - origin[2]
        )]
	tag_bbox(catalog, source_relative)
	catalog[, inside0 := inside]
	catalog[, rn := .I]
	setkey(catalog, rn)

	if (catalog[, any(inside)]) {
		# tag Inside Source
		if (tagBySourceName) {
			if (tagPoly) {
				TDinside <- source_relative[,
				{
					tag_bbox(catalog[, inside := inside0], .(x = x, y = y))
					cbind(ID = catalog[(inside), rn], pnt.in.poly(catalog[(inside), cbind(x, y)], cbind(x, y)))
				}, by = .(area, pid)][pip == 1L, paste0(paste(area, pid, sep = ","), collapse = ";"), by = ID]	
			} else {
				TDinside <- source_relative[,
				{
					tag_bbox(catalog[, inside := inside0], .(x = x, y = y))
					cbind(ID = catalog[(inside), rn], pnt.in.poly(catalog[(inside), cbind(x, y)], cbind(x, y)))
				}, by = .(area, pid)][pip == 1L, paste0(area, collapse = ","), by = ID]					
			}
			catalog[TDinside, tagInside := V1]
		} else {
			TDinside <- source_relative[,
			{
				tag_bbox(catalog[, inside := inside0], .(x = x, y = y))
				cbind(ID = catalog[(inside), rn], pnt.in.poly(catalog[(inside), cbind(x, y)], cbind(x, y)))
			}, by = .(area, pid)][, sum(pip), by = ID]
			catalog[TDinside, tagInside := as.logical(V1)]
		}
	}
	catalog[, ":="(inside = NULL, inside0 = NULL, rn = NULL)]

	invisible(catalog)
}


