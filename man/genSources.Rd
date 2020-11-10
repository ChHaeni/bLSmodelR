\name{genSources}
\alias{genSources}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Source Geometry
}
\description{
Generate a \code{data.frame} providing source geometry data.
}
\usage{
genSources(\dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
a \code{list}, a \code{data.frame} and/or individual arguments defining source areas. Names of list entries and/or arguments will be taken as source names. (see examples below)
}

}
\details{
see examples below on the usage of \code{genSources}...
}
\value{
a \code{data.frame} of class \code{Sources}
}
\author{
Christoph Haeni
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{Class-Sources}}, \code{\link{genInputList}}, \code{\link{genInterval}}, \code{\link{genSensors}}, \code{\link{genModel}}, \code{\link{genTolerances}}, \code{\link{runbLS}}, \code{\link{bLSmodelR-package}}.
}
\examples{
\dontrun{
## generate 2 circular source areas:
Sources1 <- genSources(list(
  "Circle1" = list("c", M = c(0,0), R = 10),
  "Circle2" = list("c", M = c(0,30), R = 20, N = 50))
  )
head(Sources1)

## generate a rectangular and a polygon source area:
Sources2 <- genSources(
  RectSource = list("r", x1 = -10, x2 = 10, y1 = 50, y2 = 70),
  Polygon = list(x = c(-40, -30, -20, -30), y = c(30, 20, 30, 40))
  )
Sources2

# join Sources:
Sources <- join(Sources1,Sources2)
head(Sources)
## look at site map:
siteMap(Sources,PolyArgs=list(col="darkgreen"),SourceTextArgs=list(cex=1,col="navyblue"))
addScaleBar()


## generate a source consisting of 3 circular polygon areas:
SouList <- lapply(list(c(-20,5),c(-20,-5),c(-11.3,0)),function(x)list("c",M=x,R=1,N=50))
names(SouList) <- rep("ArtSource",3)
Sources3 <- genSources(SouList)
head(Sources3)

siteMap(Sources3,PolyArgs=list(col="darkgreen")
	,SourceTextArgs=list(cex=1,col="navyblue"),xlim=c(-5,15),ylim=c(-10,10))
addScaleBar()


## generate a rectangular source area from a data frame:
SouDF <- data.frame("Rect",c(-10,-10,10,10),c(50,70,70,50),1)
Sources4 <- genSources(SouDF)
Sources4

plot(join(Sources,Sources3,Sources4))

}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
