\name{genTolerances}
\alias{genTolerances}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation Tolerances
}
\description{
Generate a \code{data frame} defining calculation tolerances.
}
\usage{
genTolerances(Tol.Z = 2.5, Tol.L = 10, Tol.Zo = 5, Tol.sUu = 10, Tol.sVu = 10, 
  Tol.sWu = 5, Tol.Zero = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tol.Z}{
Tolerance (in \%) of sensor height when checking existing touchdown catalogs.
}
  \item{Tol.L}{
Tolerance (in \%) of Obukhov-Length when checking existing touchdown catalogs.
}
  \item{Tol.Zo}{
Tolerance (in \%) of roughness height when checking existing touchdown catalogs.
}
  \item{Tol.sUu}{
Tolerance (in \%) of sigmaU/ustar when checking existing touchdown catalogs.
}
  \item{Tol.sVu}{
Tolerance (in \%) of sigmaV/ustar when checking existing touchdown catalogs.
}
  \item{Tol.sWu}{
Tolerance (in \%) of sigmaW/ustar (at sensor height) when checking existing touchdown catalogs.
}
  \item{Tol.Zero}{
logical. If \code{TRUE}, all touchdown catalog selection tolerances will be set to zero.
}

}
% \details{
% some details later...
% }
% \value{
% kommt noch...
% }
% \references{

% }
\author{
Christoph Haeni
}
% \note{
% some notes later...
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{genInputList}}, \code{\link{genInterval}}, \code{\link{genSensors}}, \code{\link{genSources}}, \code{\link{model_options}}, \code{\link{runbLS}}, \code{\link{bLSmodelR-package}}.
}
\examples{
\dontrun{
## default:
genTolerances()

## Zero tolerances:
genTolerances(Tol.Zero=TRUE)
}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
