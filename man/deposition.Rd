\name{deposition}
\alias{deposition}

\title{
Dry Deposition Post-Processing 
}
\description{
Run the \sQuote{bLSmodelR} dry deposition post-processing function.
}
\usage{
deposition(x, vDep, rn = NULL, Sensor = NULL, Source = NULL, vDepSpatial = NULL, ncores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A \code{data.frame} (or \code{data.table}) of class \code{bLSresult} output by \code{\link{runbLS}} containing the results of a model run.
}
  \item{vDep}{
either a numeric vector providing the model surface deposition velocity (at height d + z0) for the rows specified by argument \code{index} or a character string specifying the column that contains the \code{vDep} values.
}
  \item{rn}{
integer. Subsetting by column \sQuote{rn}. Only matching rows will be post-processed.
}
  \item{Sensor}{
character. Subsetting by column \sQuote{Sensor}. Only matching rows will be post-processed.
}
  \item{Source}{
character. Subsetting by column \sQuote{Source}. Only matching rows will be post-processed.
}
  \item{vDepSpatial}{
a \code{list} with two entries (this argument is still experimental!). First \code{list} entry is a named \code{list} or named \code{vector} that defines area names with corresponding deposition velocities. Second \code{list} entry is an object of class \code{Sources} (see \code{\link{genSources}}) that defines the geometry of those areas (see example below). Note: Areas should not overlap each other!
}
  \item{ncores}{
integer. Number of cores used for parallel processing. 
}
}
% \details{
% %The results of the model run can be saved as csv file, by... 
% }
\value{
A \code{data.frame} (or \code{data.table}) of class \code{bLSresult} and \code{deposition} containing the results of the model run and the deposition post-processing.
The following \code{attribute} is additionally attached to the existing \code{bLSresult} \code{attributes}:
\describe{
\item{vDep}{vector of surface deposition velocities}
}

}

\author{
Christoph Haeni
}
% \note{
% %%
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~
% \references{
% Flesch, T. K., J. D. Wilson, et al. (1995). \dQuote{Backward-time Lagrangian stochastic dispersion models and their application to estimate gaseous emissions.} Journal of Applied Meteorology 34(6): 1320-1332.
% }

\seealso{
\code{\link{bLSmodelR-package}}, \code{\link{genInputList}}, \code{\link{coreModel}}.
}
\examples{
\dontrun{
# Example Run with default values:
Sensor <- genSensors(
  PointSensor1=list(x=0,y=0,z=2)
  ,PointSensor2=list(x=10,y=0,z=2)
  )
Sources <- genSources("Circle 1" = list("c", M = c(0,50), R = 10))
Ints <- genInterval(MaxFetch=70)
InputList <- genInputList(Sensor,Ints,Sources)

## run model
Run <- runbLS(InputList,Cat.Path=getwd())

## dry deposition post-processing:
RunDep <- deposition(Run,vDep=0.035,Sensor='PointSensor1')

## spatially inhomogeneous deposition velocity (still experimental!)
DepAreas <- genSources(
	"Circle 2" = list("c", M = c(0,30), R = 10)
	,"Circle 3" = list("c", M = c(0,10), R = 10)	
)
DepVel <- list("Circle 2" = 0.005,"Circle 3" = 0.01)
plot(Sensor,join(Sources,DepAreas),SourceTextArgs=list(labels=paste0('vDep = ',c(0,0.005,0.01))))
RunDep2 <- deposition(Run,vDep=0.035,Sensor='PointSensor1',vDepSpatial=list(DepVel,DepAreas))


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
