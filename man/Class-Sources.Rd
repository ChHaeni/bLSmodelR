\name{Class-Sources}
\alias{Class-Sources}
\alias{Sources}
\alias{print.Sources}
\alias{head.Sources}
\alias{tail.Sources}
\alias{join.Sources}
\alias{plot.Sources}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S3 methods for class \code{Sources}
}
\description{
S3 methods for class \code{Sources}:
}

\usage{
\method{print}{Sources}(x, Nrows = 30, \dots)
\method{head}{Sources}(x, \dots)
\method{tail}{Sources}(x, \dots)

\method{plot}{Sources}(x, \dots)

\method{join}{Sources}(\dots)

}

\arguments{
  \item{x}{
a \code{data.frame} of class \code{Sources}.
}	
  \item{Nrows}{
number of maximum rows to print without changing to compact display.
}	
  \item{\ldots}{
further arguments to be passed to other methods. For \code{join.Sources} see the details section below.
}	
}

\details{
For \code{join.Sources} no further arguments other than different \code{data.frame}s of class \code{Sources} can be supplied. 
}


\author{
Christoph Haeni
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bLSmodelR-package}}, \code{\link{runbLS}}, \code{\link{genSources}}
}
