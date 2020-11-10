\name{Class-TDcat}
\alias{Class-TDcat}
\alias{TDcat}
\alias{print.TDcat}
\alias{plot.TDcat}
\alias{points.TDcat}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{S3 methods for class \code{TDcat}}

\description{
S3 methods for class \code{TDcat}:
}

\usage{
\method{print}{TDcat}(x, \dots)

\method{plot}{TDcat}(x, asp = 1, panel.first = {
    grid()
    abline(h = 0, col = "darkgrey")
    abline(v = 0, col = "darkgrey")
},\dots)
\method{points}{TDcat}(x, \dots)

}

\arguments{
  \item{x}{
a \code{data.table} of class \code{TDcat}.
}	
  \item{asp}{
y/x aspect ratio. see \code{\link{plot.window}}.
}	
  \item{panel.first}{
plotting expression prior to the main plotting. see \code{\link{plot.default}}.
}	
  \item{\ldots}{
further arguments to be passed to other methods.
}	
}


\author{
Christoph Haeni
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bLSmodelR-package}}, \code{\link{runbLS}}, \code{\link{coreModel}}
}
