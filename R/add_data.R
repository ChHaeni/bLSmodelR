add_data <- function(x, add, on = intersect(names(x), names(add))){
    temp <- copy(x)
    setDT(temp)
    out <- merge(temp, add, by = on, all.x = TRUE)
    # add attributes
    add_attr <- names(attributes(x)) %w/o% c(".internal.selfref", "row.names", "names")
    for(att in add_attr){
        setattr(out, att, attr(x, att))
    }
    out
}
