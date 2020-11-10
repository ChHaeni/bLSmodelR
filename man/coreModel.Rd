\name{coreModel}
\alias{coreModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate trajectories in backward mode.
}
\description{
This is the core function to calculate the (air parcel) trajectories in a backward mode (i.e. backward-in-time).
}
\usage{
coreModel(u, v, w, zSens, ustar, L, Zo, bw, sigmaUustar, sigmaVustar, 
  kv, C0, alpha, MaxFetch, method = "Euler")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{u}{
  numeric vector. Instantaneous velocity of the individual air parcels at time 0 in alongwind direction (given in m/s).
  }
    \item{v}{
  instantaneous velocity of air parcel at time 0 in crosswind direction (given in m/s).
  }
    \item{w}{
  instantaneous velocity of air parcel at time 0 in vertical direction (given in m/s).
  }
    \item{zSens}{
  sensor height in m above ground level.
  }
    \item{ustar}{
  friction velocity in m/s.
  }
    \item{L}{
  Obukhov length in m.
  }
    \item{Zo}{
  roughness length in m.
  }
    \item{bw}{
  parameter defining the vertical profile of the variance in the vertical velocity component sigmaW. 
  }
    \item{sigmaUustar}{
  variance of the alongwind velocity devided by \code{ustar}.
  }
    \item{sigmaVustar}{
  variance of the crosswind velocity devided by \code{ustar}.
  }
    \item{kv}{
  von Ka'rma'n constant (usually = 0.4).
  }
    \item{C0}{
  Kolmogorov constant.
  }
    \item{alpha}{
  Fraction of the velocity decorrelation time scale as given in Flesch et al., 2004 to choose the time increment.
  }
    \item{MaxFetch}{
  maximum tracking distance of trajectories (in m).
  }
    \item{method}{
  the method of numerical discretization to be used. The default \dQuote{Euler} corresponds to \cite{Flesch et al. (1995)}. \dQuote{RK4} corresponds to the classical (4th order) Runge-Kutta method.
  }
}
% \details{
% einige details...
% The profiles of the velocity statistics are taken from MOST, as described in ...
% }
\value{
A list with following items:
\item{Traj_IDOut}{ID of the corresponding trajectory. This ID is identical to the indexing position of the initial velocities \code{u}, \code{v} and \code{w} (see calculation of w'C'/E in example below).}
\item{TimeOut}{Time until touchdown in seconds.}
\item{xOut}{x position of touchdown.}
\item{yOut}{y position of touchdown.}
\item{wTDOut}{Touchdown velocity.}
}
\references{
Flesch, T. K., J. D. Wilson, et al. (1995). \dQuote{Backward-time Lagrangian stochastic dispersion models and their application to estimate gaseous emissions.} Journal of Applied Meteorology 34(6): 1320-1332.

Flesch, T. K., J. D. Wilson, et al. (2004). \dQuote{Deducing ground-to-air emissions from observed trace gas concentrations: A field trial.} Journal of Applied Meteorology 43(3): 487-502.
}
\author{
Christoph Haeni
}
% \note{
% some notes. coreModelWrapper...
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{runbLS}}, \code{\link{bLSmodelR-package}}.
}
\examples{
\dontrun{

## set up model parameters
zSigmaWu <- 2
zSens <- 2
Ustar <- 0.25
L <- -600
Zo <- 0.001
SigmaUu <- 4.5
SigmaVu <- 4
SigmaWu <- 1.7
kv <- 0.4
alpha <- 0.02
MaxFetch <- 500
A <- 0.5

bw <- calcbw(SigmaWu, zSigmaWu/L)
SigmaWm <- calcsigmaW(Ustar, zSens/L, bw)
U <- calcU(Ustar, Zo, L, zSens, kv)
C0 <- calcC0(bw, kv)

## initial velocities
N0 <- 1000
set.seed(1234)
TDcat <- initializeCatalog(N0, Ustar, SigmaUu, SigmaVu, bw,
  zSens, L, Zo, A, alpha, MaxFetch, kv)
      
uvw <- uvw0(TDcat)

## run trajectory calculation
TDList <- coreModel(uvw[,"u0"], uvw[,"v0"], uvw[,"w0"], zSens, Ustar, 
  L, Zo, bw, SigmaUu, SigmaVu, kv, C0, alpha, MaxFetch)

## assign touchdowns to initalized Catalog:
assignCat(TDcat,TDList)

## Calculate C/E and w'C'/E ratios of entire domain with homogeneous
## emission rate:
Ci <- data.table(ID=1:N0,CE=0,u=uvw[,"u0"],v=uvw[,"v0"],w=uvw[,"w0"],key="ID")
Ci[TDcat[,sum(2/wTD),by=Traj_ID],CE:=V1]

stats <- Ci[,cbind(
  round(rbind(ce <- mean(CE),se <- sd(CE)/sqrt(N0),ce - 1.96*se,ce + 1.96*se),1),
  round(rbind(wce <- mean(w*(cs <- CE-ce)),wse <- sd(w*cs)/sqrt(N0)
  ,wce - 1.96*wse,wce + 1.96*wse),2))]

dimnames(stats) <- list(c("E(X)","SE","Lower-CI95\%","Upper-CI95\%"),c("C/E","w'C'/E"))
print(stats)



## look at Catalog:
TDcat
plot(TDcat,pch=20)


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
