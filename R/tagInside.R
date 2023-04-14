# TODO:
#   - change to consistent tagging (character instead of character/logical)
#       -> careful rattail!
#   - rename function 
#       -> rename all existing tagInside() calls
#   - better arguments/argument names
#   - use better names for Sensors class (simpleNames)
#   - use tag_inside in calcCE.R ?
#   - use methods?
tagInside <- function(tag_catalog, tag_source, 
    sensor_position = c("x-Coord (m) " = 0, "y-Coord (m) " = 0),
    tagBySourceName = tagPoly, tagPoly = FALSE) {
	tag_catalog[, tagInside := 
        if (tagBySourceName) {
            rep("", nrow(tag_catalog))
        } else {
            rep(FALSE, nrow(tag_catalog))
        }
    ]
	source_relative <- data.table(tag_source)
	setnames(source_relative, c("area", "x", "y", "pid"))
	source_relative[, ":="(
        x = x - sensor_position[, "x-Coord (m) "],
        y = y - sensor_position[, "y-Coord (m) "]
        )]
	tag_bbox(tag_catalog, source_relative)
	tag_catalog[, inside0 := inside]
	tag_catalog[, rn := .I]
	setkey(tag_catalog, rn)

	if (tag_catalog[, any(inside)]) {
		# tag Inside Source
		if (tagBySourceName) {
			if (tagPoly) {
				TDinside <- source_relative[,
				{
					tag_bbox(tag_catalog[, inside := inside0], .(x = x, y = y))
					cbind(ID = tag_catalog[(inside), rn], pnt.in.poly(tag_catalog[(inside), cbind(x, y)], cbind(x, y)))
				}, by = .(area, pid)][pip == 1L, paste0(paste(area, pid, sep = ","), collapse = ";"), by = ID]	
			} else {
				TDinside <- source_relative[,
				{
					tag_bbox(tag_catalog[, inside := inside0], .(x = x, y = y))
					cbind(ID = tag_catalog[(inside), rn], pnt.in.poly(tag_catalog[(inside), cbind(x, y)], cbind(x, y)))
				}, by = .(area, pid)][pip == 1L, paste0(area, collapse = ","), by = ID]					
			}
			tag_catalog[TDinside, tagInside := V1]
		} else {
			TDinside <- source_relative[,
			{
				tag_bbox(tag_catalog[, inside := inside0], .(x = x, y = y))
				cbind(ID = tag_catalog[(inside), rn], pnt.in.poly(tag_catalog[(inside), cbind(x, y)], cbind(x, y)))
			}, by = .(area, pid)][, sum(pip), by = ID]
			tag_catalog[TDinside, tagInside := as.logical(V1)]
		}
	}
	tag_catalog[, ":="(inside = NULL, inside0 = NULL, rn = NULL)]

	invisible(tag_catalog)
}


