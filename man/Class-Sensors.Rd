\name{Class-Sensors}
\alias{Class-Sensors}
\alias{Sensors}
\alias{print.Sensors}
\alias{head.Sensors}
\alias{tail.Sensors}
\alias{join.Sensors}
\alias{plot.Sensors}



%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S3 methods for class \code{Sensors}
}
\description{
S3 methods for class \code{Sensors}:
}
 	

\usage{
\method{print}{Sensors}(x, Nrows = 30, \dots)
\method{head}{Sensors}(x, \dots)
\method{tail}{Sensors}(x, \dots)

\method{plot}{Sensors}(x, \dots)

\method{join}{Sensors}(\dots)

}

\arguments{
  \item{x}{
a \code{data.frame} of class \code{Sensors}.
}	
  \item{Nrows}{
number of maximum rows to print without changing to compact display.
}	
  \item{\ldots}{
further arguments to be passed to other methods. For \code{join.Sensors} see the details section below.
}	
}

\details{
For \code{join.Sensors} no further arguments other than different \code{data.frame}s of class \code{Sensor} can be supplied. 
}


\author{
Christoph Haeni
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bLSmodelR-package}}, \code{\link{runbLS}}, \code{\link{genSensors}}
}
