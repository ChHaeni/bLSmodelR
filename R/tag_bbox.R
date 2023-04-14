#' tag touchdowns inside bbox
#'
#' tag touchdowns inside a bbox around arbitrarily defined x/y values
#'
#' @param catalog A bLSmodelR TD catalog (class TDcat).
#' @param tag_polygon Values of x and y coordinates (intended to be a polygon for later
#'  checking with tag_inside.
tag_bbox <- function(catalog, tag_polygon) {
	if (!("inside" %in% names(catalog))) catalog[, inside := TRUE]
	UseMethod("tag_bbox", tag_polygon)
}
tag_bbox.matrix <- function(catalog, tag_polygon) {
	catalog[(inside), inside := (
		x <= max(tag_polygon[, 1]) & 
		x >= min(tag_polygon[, 1]) & 
		y <= max(tag_polygon[, 2]) & 
		y >= min(tag_polygon[, 2]))
	]
}
tag_bbox.data.table <- function(catalog, tag_polygon) {
	catalog[(inside), inside := (
		x <= tag_polygon[, max(x)] & 
		x >= tag_polygon[, min(x)] & 
		y <= tag_polygon[, max(y)] & 
		y >= tag_polygon[, min(y)])
	]
}
tag_bbox.list <- function(catalog, tag_polygon) {
	catalog[(inside), inside := (
		x <= max(tag_polygon[[1]]) & 
		x >= min(tag_polygon[[1]]) & 
		y <= max(tag_polygon[[2]]) & 
		y >= min(tag_polygon[[2]]))
	]
}
