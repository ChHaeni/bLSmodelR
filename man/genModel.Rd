\name{genModel}
\alias{genModel}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model Parameters
}
\description{
Generate a \code{data.frame} defining (global) model specific data.
}
\usage{
genModel(kv = 0.4, A = 0.5, alpha = 0.02, wTDcutoff = 1E-4, TDwrite = TRUE,
  overwriteTD = TRUE, TDread = TRUE, TDonly = FALSE, ncores = 1)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kv}{
von-Karman constant.
}
  \item{A}{
scaling constant. (see \cite{Flesch et al. (2004)})
}
  \item{alpha}{
time step fraction of the velocity decorrelation time scale. note: alpha > 10^-3
}
  \item{wTDcutoff}{
cutoff velocity used as lower bound of calculated touchdown velocities.
}
  \item{TDwrite}{
logical. Should touchdown catalogs be saved?
}
  \item{overwriteTD}{
logical. Should existing touchdowns catalogs be overwritten?
}
  \item{TDread}{
logical. Check for existing touchdown catalogs?
}
  \item{TDonly}{
logical. Skip calculation of C/E etc., only simulate touchdown catalogs? Default is \code{FALSE}.
}
  \item{ncores}{
integer. Number of cores used for parallel processing. 
}
}
% \details{
% some details later...
% }
\value{
a \code{list} with the above defined model specifications. These model parameters are used \sQuote{globally}, i.e. for all intervals, in the specified model run.
}
\references{
Flesch, T. K., J. D. Wilson, et al. (2004). \dQuote{Deducing ground-to-air emissions from observed trace gas concentrations: A field trial.} Journal of Applied Meteorology 43(3): 487-502.
}
\author{
Christoph Haeni
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{genInputList}}, \code{\link{genInterval}}, \code{\link{genSensors}}, \code{\link{genSources}}, \code{\link{genTolerances}}, \code{\link{runbLS}}, \code{\link{bLSmodelR-package}}.
}
\examples{
\dontrun{
## defaults:
genModel()
}
}


% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
