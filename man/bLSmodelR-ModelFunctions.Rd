\name{ModelFunctions}
\alias{ModelFunctions}
\alias{calcU}
\alias{calcC0}
\alias{calcepsilon}
\alias{calcsigmaW}
\alias{calcTL}
\alias{calcbw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model functions
}
\description{
Functions to calculate some desciptive parameters outside of the model.
}
\usage{
calcU(ustar, Zo, L, z, kv = 0.4)
calcepsilon(ustar, L, bw, z, kv = 0.4)
calcsigmaW(ustar, zL, bw)
calcC0(bw, kv = 0.4, A = 0.5)
calcTL(sigmaW, C0, epsilon)
calcbw(sigmaWustar, zL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ustar}{
friction velocity in m/s.
}
  \item{Zo}{
roughness length in meters.
}
  \item{L}{
Obukhov-Length in meters.
}
  \item{z}{
(aerodynamic) height in meters.
}
  \item{zL}{
\code{z/L}.
}
  \item{bw}{
scaled standard deviation of the vertical wind speed component \bold{w} at \code{z = 0}.
}
  \item{kv}{
von-Karman constant.
}
  \item{A}{
scaling constant. (see \cite{Flesch et al. (2004)})
}
  \item{sigmaW}{
standard deviation of the vertical wind speed component \bold{w}.
}
  \item{C0}{
Kolmogorov's (universal) constant
}
  \item{epsilon}{
dissipation rate of turbulent kinetic energy in m^{2}/s^{3}.
}
  \item{sigmaWustar}{
scaled standard deviation of the vertical wind speed component \bold{w}: sigmaW/ustar.
}

}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
% \value{
% %%  ~Describe the value returned
% %%  If it is a LIST, use
% %%  \item{comp1 }{Description of 'comp1'}
% %%  \item{comp2 }{Description of 'comp2'}
% %% ...
% }
\references{
Flesch, T. K., J. D. Wilson, et al. (2004). \dQuote{Deducing ground-to-air emissions from observed trace gas concentrations: A field trial.} Journal of Applied Meteorology 43(3): 487-502.
}
\author{
Christoph Haeni
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
example of \code{\link{coreModel}}, \code{\link{bLSmodelR-package}}.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
