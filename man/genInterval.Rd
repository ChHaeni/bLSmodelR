\name{genInterval}
\alias{genInterval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interval Data
}
\description{
Generate a \code{data.frame} providing data on an interval basis (i.e. processed Sonic-Anemometer data etc.).
}
\usage{
genInterval(Data = NULL, Ustar = 0.25, L = -2000, Zo = 0.01, sUu = 2.5, sVu = 2, 
  sWu = 1.25, z_sWu = 2, WD = 0, d = 0, N0=5E4, MaxFetch=500, SensorNames = "", 
  SourceNames = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
a \code{data.frame}. Columns that match argument names below will be taken as input to the correspondent argument. Columns that don't match will be appended as additional data.
}
  \item{Ustar}{
friction velocity in m/s.
}
  \item{L}{
Obukhov-Length in meters.
}
  \item{Zo}{
roughness length in meters.
}
  \item{sUu}{
scaled standard deviation of the wind speed component \bold{u} (i.e. alongwind): sigmaU/ustar.
}
  \item{sVu}{
scaled standard deviation of the wind speed component \bold{v} (i.e. crosswind): sigmaV/ustar.
}
  \item{sWu}{
scaled standard deviation of the vertical wind speed component \bold{w}: sigmaW/ustar.
}
  \item{z_sWu}{
measurement height of the wind statistics sigmaW/ustar (\code{sWu}) in meters. Only relevant for unstable conditions (\code{L}<0), where sigmaW/ustar is assumed to vary with height. Usually, this is the Sonic-Anemometer measurement height above ground.
}
  \item{WD}{
average wind direction in degrees deviating from North. I.e. wind blowing from 0/360 = north, 90 = east, 180 = south, 270 = west.
}
  \item{d}{
displacement height in meters.
}
  \item{N0}{
integer. Number of trajectories in model runs.
}
  \item{MaxFetch}{
numeric. Maximum upwind fetch (alongwind tracking distance). Negative numbers will try to optimize this value (see details below).
}
  \item{SensorNames}{
Names of the Sensors to calculate TDs, separated by \sQuote{,}.
}
  \item{SourceNames}{
Names of the Sources to include, separated by \sQuote{,}.
}

}
\details{
Providing arguments directly will overwrite corresponding input from argument Data (see examples below).

Providing a negative \code{MaxFetch} argument will calculate the maximum distance between the sensor and the source and add the absolute value of the provided \code{MaxFetch} argument to this distance (see example below). 
}
\value{
a \code{data.frame} of class \code{Sonic}.
}
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
\code{\link{genInputList}}, \code{\link{genSensors}}, \code{\link{genSources}}, \code{\link{model_options}}, \code{\link{genTolerances}}, \code{\link{runbLS}}, \code{\link{bLSmodelR-package}}.
}
\examples{
\dontrun{
## default:
print(x <- genInterval())
class(x)

## providing Data argument:
dat <- data.frame(Ustar=0.2,L=c(-30,-300),add1='test',ustar=0.1)
genInterval(dat)

## partially overwriting Data argument:
genInterval(dat,Ustar=0.3)

## optimizing MaxFetch
Sensor <- genSensors(PS=list(x=0,y=0,z=2))
Source <- genSources(
  Src = list('c',M=c(0,50),R=20)
  ,Src = list(x = c(25,35,35,25), y = c(20,20,50,50))
  )
Ints <- genInterval(WD=c(0,22.5,45,67.5),MaxFetch=-20)
InList <- genInputList(Ints,Source,Sensor)

IntExt <- prepareIntervals(InList,getwd())
IntExt[,MaxFetch]

## checking graphically:
plot(InList)
IntExt[,{
  x0 <- sin(WD/180*pi)
  y0 <- cos(WD/180*pi)
  x <- x0*(MaxFetch - 20)
  y <- y0*(MaxFetch - 20)
  for(i in 1:4){
    lines(c(0,x[i]),c(0,y[i]),lty=i)
    lines(x[i]+c(-100,100)*y0[i],y[i]-c(-100,100)*x0[i],lwd=2,col='lightgrey',lty=i)
  }}]

}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

