#' tag touchdowns inside bbox
#'
#' tag touchdowns inside a bbox around arbitrarily defined x/y values. Adds a column
#'  "bbox_inside" to the TDcat object
#'
#' @param catalog A TDcat object, i.e., a bLSmodelR TD catalog of class TDcat.
#' @param tag_polygon Values of x and y coordinates (intended to be a polygon for later
#'  checking with tag_inside.
#' @return The provided TDcat object with a column named "bbox_inside" indicating if
#'  touchdowns are inside the bbox or not (TRUE/FALSE).
tag_bbox <- function(catalog, tag_polygon) {
	if (!("bbox_inside" %in% names(catalog))) catalog[, bbox_inside := TRUE]
	UseMethod("tag_bbox", tag_polygon)
}
tag_bbox.matrix <- function(catalog, tag_polygon) {
	catalog[(bbox_inside), bbox_inside := (
		x <= max(tag_polygon[, 1]) & 
		x >= min(tag_polygon[, 1]) & 
		y <= max(tag_polygon[, 2]) & 
		y >= min(tag_polygon[, 2]))
	]
}
tag_bbox.data.table <- function(catalog, tag_polygon) {
	catalog[(bbox_inside), bbox_inside := (
		x <= tag_polygon[, max(x)] & 
		x >= tag_polygon[, min(x)] & 
		y <= tag_polygon[, max(y)] & 
		y >= tag_polygon[, min(y)])
	]
}
tag_bbox.list <- function(catalog, tag_polygon) {
	catalog[(bbox_inside), bbox_inside := (
		x <= max(tag_polygon[[1]]) & 
		x >= min(tag_polygon[[1]]) & 
		y <= max(tag_polygon[[2]]) & 
		y >= min(tag_polygon[[2]]))
	]
}
