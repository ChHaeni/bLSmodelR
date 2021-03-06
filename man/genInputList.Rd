\name{genInputList}
\alias{genInputList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model Input
}
\description{
Create a \code{list} gathering the required model input.
}
\usage{
genInputList(..., Tol.Zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\ldots}{
optional arguments containing pre-defined model input: 
\describe{
    \item{\code{Class-Interval}}{a \code{data.frame} as generated by \code{\link{genInterval}} containing the Interval input.}
    \item{\code{\link{Class-Sensors}}}{a \code{data.frame} as generated by \code{\link{genSensors}} containing the Sensors input.}
    \item{\code{\link{Class-Sources}}}{a \code{data.frame} as generated by \code{\link{genSources}} containing the Sources input.}
    \item{\code{Class-Model}}{a \code{list} as generated by \code{\link{genModel}} containing the Model parameter input.}
    \item{\code{Class-Tolerances}}{a \code{data.frame} as generated by \code{\link{genTolerances}} containing the Tolerances input.}
}
}
  \item{Tol.Zero}{
logical. If \code{TRUE}, all touchdown catalog selection tolerances will be set to zero.
}

}
% \details{
% some details later...
% }
\value{
a \code{list} containing all necessary model input.
}
% \references{
% Flesch, T. K., J. D. Wilson, et al. (2004). \dQuote{Deducing ground-to-air emissions from observed trace gas concentrations: A field trial.} Journal of Applied Meteorology 43(3): 487-502.
% }
\author{
Christoph Haeni
}
% \note{
% some notes later...
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{genInterval}}, \code{\link{genSensors}}, \code{\link{genSources}}, \code{\link{genModel}}, \code{\link{genTolerances}}, \code{\link{runbLS}}, \code{\link{bLSmodelR-package}}.
}
\examples{
\dontrun{
## generate source areas:
Sources <- genSources("Circle 1" = list("c", M = c(0,0), R = 10),
  "Circle 2" = list("c", M = c(0,30), R = 20, N = 50))

## generate sensors:
Sensors <- genSensors(LineSensor1=list(x=-30,y=c(-10,50),z=1.25,n=40)
  ,LineSensor2=list(x=30,y=c(-10,50),z=c(1.05,1.25),d=0.5))

## look at site map:
siteMap(Sensors,Sources,PolyArgs=list(col="darkgreen")
  ,PtArgs=list(pch=3,cex=0.8,lwd=2),SourceTextArgs=list(cex=1,col="navyblue")
  ,SensorTextArgs=list(cex=0.7,pos=3),xlim=c(-60,60))
addScaleBar()
addWindrose(WD = 270)

# generate interval data:
Int <- genInterval(WD=270,N0=1E4,SensorNames="LineSensor2",SourceNames="Circle 1,Circle 2")

## generate List of model input:
Model <- genModel(TDwrite=FALSE)
InputList <- genInputList(Int, Sources, Sensors, Model)

## run bLSmodelR (takes 2 to 5 mins):
Run <- runbLS(InputList,Cat.Path=getwd())

# switch names to "informative" and extract specific results (Source column is sorted accordingly)
switchNames(Run,simple=FALSE)
extractResult(Run,list(Source=c("Circle 2","Circle 1")),keep=c("Sensor","Source","C/E","C/E SE"))

## clean up TD catalogs
cleanTemporary()
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
