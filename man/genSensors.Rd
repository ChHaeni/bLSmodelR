\name{genSensors}
\alias{genSensors}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sensor Geometry
}
\description{
Generate a \code{data.frame} providing sensor geometry data.
}
\usage{
genSensors(\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
a \code{data.frame} and/or individual arguments defining point and/or line sensors. Names of arguments will be taken as sensor names. (see examples below)
}

}
\details{
see examples below on the usage of \code{genSensors}...
}
\value{
a \code{data.frame} of class \code{Sensors}
}

\author{
Christoph Haeni
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Class-Sensors}}, \code{\link{genInputList}}, \code{\link{genInterval}}, \code{\link{genSources}}, \code{\link{genModel}}, \code{\link{genTolerances}}, \code{\link{runbLS}}, \code{\link{bLSmodelR-package}}.
}
\examples{
\dontrun{
## generate a point and two line sensors (way1):
Sensors1 <- genSensors(PointSensor=list(x=0,y=0,z=1.5)
  ,LineSensor1=list(x=c(-10,10),y=0,z=1.25,n=40)
  ,LineSensor2=list(x=0,y=c(-20,10),z=c(1.05,1.25),d=0.5))
Sensors1

## generate a point and two line sensors (way2):
SensorDF2 <- data.frame(c("PointSensor2","LS1b_1","LS1b_2","LS2b_1","LS2b_2"),c(0,-10,10,0,0)
  ,c(0,0,0,-20,10),c(1.5,1.25,1.25,1.05,1.25),c("","LineSensor1b","LineSensor1b","LineSensor2b",
  "LineSensor2b"),c(1,40,0,NA,NA),c(0,0,0,0.5,0))
Sensors2 <- genSensors(SensorDF2)
Sensors2

## join Sensors data frames:
Sensors4 <- join(Sensors1,Sensors2)
Sensors4

siteMap(Sensors4,PtArgs=list(pch=4,cex=0.8,lwd=2,col="darkred")
  ,SensorTextArgs=list(cex=0.7,pos=3),LSArgs = list(col="navyblue"))
addScaleBar()

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
