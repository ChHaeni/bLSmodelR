# TODO:
# -> join if sensor name already exists? as argument?
# -> rename_source

#' Rename sensors
#'
#' @param x object where sensors should be renamed
#' @param old character. Old sensor names.
#' @param new character. New sensor names.
#' @param \dots further arguments which are passed to methods. Ignored.
#' @return the original object with modified sensor names
#' @export
#'
rename_sensor <- function(x, old, new) {
    # check if old is character
    if (length(old) < 1 || !is.character(old)) {
        stop('argument "old" is not valid')
    }
    # check if new is character
    if (length(new) < 1 || !is.character(new)) {
        stop('argument "new" is not valid')
    }
    # check equal lengths old = new
    if (length(old) != length(new)) {
        stop('arguments "old" and "new" must be of equal length')
    }
    # check if valid names
    # TODO: allow below
	if(any(bad_new <- grepl("[.][0-9]*$", new))){
        stop("Sensor Names are not allowed to end with .[0-9]*!\n\t -> Sensor Names: ",
            paste(paste0('"', new[bad_new], '"'), collapse = ", "),
            " are not valid!\n")
	}
    # convert if version <4.2+
    x <- convert(x)
    # call method -> default => warning
    UseMethod('rename_sensor')
}

rename_sensor.bLSresult <- function(x, old, new, ...) {
    # leave original x as is
    x <- copy(x)
    # prepare key
    names(new) <- old
    # fix ModelInput & check old
    attr(x, 'ModelInput') <- rename_sensor.InputList(
        attr(x, 'ModelInput'), old, new, ...)
    # fix results
    if (switch_df <- !is.data.table(x)) setDT(x)
    x[Sensor %chin% old, Sensor := new[Sensor]]
    if (switch_df) setDF(x)
    # fix CalcSteps
    attr(x, 'CalcSteps')[Sensor %chin% old, c('Sensor', 'Calc.Sensor') := {
        new_cs <- Calc.Sensor
        for (o in old) {
            new_cs[Sensor == o] <- gsub(o, new[o], new_cs[Sensor == o], fixed = TRUE)
        }
        .(
            new[Sensor],
            new_cs
        )
    }]
    # fix Catalogs
    attr(x, 'Catalogs')[Sensor %chin% old, c('Sensor', 'PointSensor') := {
        new_ps <- PointSensor
        for (o in old) {
            new_ps[Sensor == o] <- gsub(o, new[o], new_ps[Sensor == o], fixed = TRUE)
        }
        .(
            new[Sensor],
            new_ps
        )
    }]
    # return
    x
}

rename_sensor.InputList <- function(x, old, new, ...) {
    # fix Sensors (and check if old exists)
    x[['Sensors']] <- rename_sensor.Sensors(x[['Sensors']],
        old, new, ...)
    # fix Interval
    old_pattern <- paste0('\\b', old, '\\b')
    for (i in seq_along(old_pattern)) {
        x[['Interval']][, 'Sensor Names (sep = ",")'] <-
            gsub(old_pattern[i], new[i], 
                x[['Interval']][, 'Sensor Names (sep = ",")'])
    }
    # return
    x
}

rename_sensor.Sensors <- function(x, old, new, ...) {
    # check if old exists in x
    if (!all(old_exists <- old %in% x[[1]])) {
        stop('Old sensor name(s): ', 
            paste(paste0('"', old[!old_exists],
                    '"'), collapse = ', '),
            ' do not exist')
    }
    # check if new already exists
    if (any(new_exists <- new %in% x[[1]])) {
        stop('New sensor name(s): ', 
            paste(paste0('"', new[new_exists],
                    '"'), collapse = ', '),
            ' does already exist. This may be allowed in a future package version.')
    }
    # prepare key
    names(new) <- old
    # replace
    x[[1]][x[[1]] %in% old] <- new[x[[1]][x[[1]] %in% old]]
    # return
    x
}
