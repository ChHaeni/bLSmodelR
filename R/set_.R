setDT.bLSresult <- function(x, keep.rownames = FALSE, key = NULL, check.names = FALSE){
	if(!is.data.table(x)){
        name = substitute(x)
        home <- function(x, env) {
            if (identical(env, emptyenv())) 
                stop("Can not find symbol ", cname, call. = FALSE)
            else if (exists(x, env, inherits = FALSE)) 
                env
            else home(x, parent.env(env))
        }

        atts <- attributes(x)
        if(check.names)x <- switchNames(x)

        cname = as.character(name)
        envir = home(cname, parent.frame())
        if (bindingIsLocked(cname, envir)) {
            stop("Can not convert '", cname, "' to data.table by reference because binding is locked. It is very likely that '", 
                cname, "' resides within a package (or an environment) that is locked to prevent modifying its variable bindings. Try copying the object to your current environment, ex: var <- copy(var) and then using setDT again.")
        }

        rn = if (!identical(keep.rownames, FALSE)) 
            rownames(x)
        else NULL
        setattr(x, "row.names", .set_row_names(nrow(x)))
        if (check.names) 
            setattr(x, "names", make.names(names(x), unique = TRUE))
    	for(i in c("CalcSteps","CatPath","Catalogs","ModelInput","sessionInfo","ModelRunTime"))setattr(x,i,atts[[i]]) 
        setattr(x, "class", c("bLSresult","data.table","data.frame"))
        alloc.col(x)
        if (!is.null(rn)) {
            nm = c(if (is.character(keep.rownames)) keep.rownames[1L] else "rn", 
                names(x))
            x[, `:=`((nm[1L]), rn)]
            setcolorder(x, nm)
        }

        if (!is.null(key))setkeyv(x, key)

        name = as.character(name)
        assign(name, x, parent.frame(), inherits = TRUE)
    }
    invisible(x)

}
setDF.bLSresult <- function(x, rownames = NULL){
	if(is.data.table(x)){
        atts <- attributes(x)
     	if (any(duplicated(rownames))) 
            stop("rownames contains duplicates")
        
        if (is.null(rownames)) {
            rn <- .set_row_names(nrow(x))
        }
        else {
            if (length(rownames) != nrow(x)) 
                stop("rownames incorrect length; expected ", 
                  nrow(x), " names, got ", length(rownames))
            rn <- rownames
        }
    	for(i in c("CalcSteps","CatPath","Catalogs","ModelInput","sessionInfo","ModelRunTime"))setattr(x,i,atts[[i]]) 
        setattr(x, "row.names", rn)
        setattr(x, "class", c("bLSresult","data.frame"))
        setattr(x, "sorted", NULL)
        setattr(x, ".internal.selfref", NULL)
    }
	invisible(x)		

}
setDT <- function(x, keep.rownames = FALSE, key = NULL, check.names = FALSE)UseMethod("setDT",x)
setDT.default <- data.table::setDT
setDF <- function(x, rownames = NULL)UseMethod("setDF",x)
setDF.default <- data.table::setDF
