\name{bLSmodelR-package}
\alias{bLSmodelR-package}
\alias{bLSmodelR}
\docType{package}
\title{
bLSmodelR - An atmospheric dispersion model in R
}
\description{
The package \bold{bLSmodelR} provides functions to set up and run a backward Lagrangian Stochastic (bLS) model.
}

% \note{
% %
% }

\author{
Christoph Haeni
}

\section{bLS Model}{
The model is a first-order Lagrangian Stochastic model that is run in backward mode (i.e. backward in time), assuming horizontally homogeneous and vertically inhomogeneous Gaussian turbulence. It is based on the paper published by \cite{Flesch et al. (2004)}. The basis of the model calculation is a generalized Langevin equation given as \deqn{du_i = a_i(x, u)*dt + b_{ij}(x, u)*d\xi_j}{du_i = a_i(x, u)*dt + b_i_j(x, u)*d\xi_j} \deqn{dx_i = u_i*dt}{dx_i = u_i*dt ,} where \eqn{d\xi} in the noise term is a random increment from a Gaussian distribution \eqn{N(0,dt)}. The indices \eqn{i, j = 1, 2, 3} represent the alongwind (x), crosswind (y) and vertical axis (z) of the rotated coordinate system. \eqn{x_i} and \eqn{u_i} are position (x, y, z) and velocity (u, v, w) of the trajectory at time t. An ensemble of upwind trajectories is calculated starting from the sensor position \eqn{x = (0, 0, z_{meas} - d)}{x = (0, 0, z_meas - d)} and the position and velocity of trajectory touchdowns on ground (i.e. on a horizontal plane at \eqn{z = z_0 + d}{z = Zo + d}, where the modelled, average wind speed equals 0.) are recorded to calculate the sensor concentration to source strength relationship as the sum of the inverse touchdown velocities, summed over all touchdowns within the source \deqn{\frac{C}{E} = \frac{2}{N_{trj}}*\sum_{src}\frac{1}{w_{TD}}}{C/E = 2/N_traj*Sum_allTDinSource(1/w_TD) .} Source areas are defined as homogenously emitting sources, having emission rate \eqn{E}. As a consequence, the \eqn{C/E} relationship is expressed in \emph{s/m}. In analogy, the flux to source strength ratios (\eqn{\left\langle u'C'\right\rangle/E, \left\langle v'C'\right\rangle/E, \left\langle w'C'\right\rangle/E}{<u'C'>/E, <v'C'>/E, <w'C'>/E}) are calculated, mainly for completeness and for possible vertical flux footprint corrections. Trajectories touching ground are perfectly reflected and their velocities are reversed to maintain covariances past reflection. Therefore, no deposition is modelled when trajectory touchdown happens. Further, this model is not capable of modelling chemistry, and trace gases have to be assumed to be reasonably inert within the relevant travelling time. For further details on the model and the concentration to emission rate relationship, see \cite{Flesch et al. (2004)}, \cite{Flesch (1996)} and  partly \cite{Flesch et al. (1995)}.
}

\section{MOST Profiles}{
The model calculation is based on Monin-Obukhov similarity theory (MOST) profiles of the wind speed and the wind statistics. As a consequence, the wind characteristics in the model domain can be set up completely by providing the three parameters
\eqn{ustar} (the friction velocity), \eqn{L} (the Obukhov length) and \eqn{z_0}{Zo} (the roughness length). It should be kept in mind that, if MOST can not be applied reasonably, model calculations certainly will be erroneous.

\emph{Wind Speed}

The vertical profile of the ensemble average of the horizontal wind speed \eqn{U = \left\langle u\right\rangle}{U = <u>} is defined as
\deqn{U(z) = ustar/k_v*\left\{ln(z/z_0) + \Psi(z/L) - \Psi(z_0/L)\right\}}{U(z) = ustar/kv*{ln(z/Zo) + \Psi(z/L) - \Psi(Zo/L)} ,}
where \eqn{\Psi} is a stability correction function defined as \deqn{if(L >= 0), \Psi(x) = 4.8*z/L} \deqn{if(L < 0), \Psi(x) = -2*ln((1+\alpha)/2) - ln((1+\alpha^2)/2) - 2*atan(\alpha) + \pi/2 ,} where \deqn{\alpha = (1 - 16*z/L)^{1/4}}{\alpha = (1 - 16*z/L)^(1/4)}

\emph{Wind Statistics}

The variances of the different velocity components are assumed to be constant within the model domain (i.e. homogenous), with the exception of the vertical velocity component under unstable atmospheric conditions (i.e. if the Obukhov length is negative, \eqn{L < 0}), where the vertical profile function is defined as \deqn{\sigma w^2 = ustar^2*bw^2*(1 - 3*z/L)^{2/3}}{\sigmaw^2 = ustar^2*bw^2*(1 - 3*z/L)^(2/3) ,} and \emph{bw} is derived from the supplied \eqn{\sigma w/ustar}{\sigmaw/ustar} ratio at given height.

\emph{Zero-Plane Displacement / Displacement Height (d)}

In contrast to \cite{Flesch et al. (2004)}, the zero-plane displacement (\eqn{d}) is included in the model profiles, by subtracting d from any measured, geometric height to give the corresponding aerodynamic heights \eqn{z = z_{meas} - d}{z = z_meas - d}.
} 

\section{Initiaization of velocities}{
Initial velocities are calculated by using an orthogonal projection of random ... (explain further)	
}

\section{Model Input}{

All model input is supplied from script/console. See \code{\link{genInputList}} for further information.

\emph{Touchdown Catalogs}

Existing touchdown catalogs can significantly speed up model runs.

\emph{Tolerances}

If touchdown catalogs will be read, the tolerances for touchdown selection can be set.
}

\section{Model Output}{

\emph{csv File}

The results of the model run can be saved to a csv file (semicolon separated).

\emph{data.frame/data.table}

All results will be given as \code{data.frames} or \code{data.table}. See \code{\link{runbLS}} for further details.

\emph{Touchdown Catalogs}

Touchdown catalogs can be saved for future runs.

}

\section{Parallel Computing}{
This package provides the option to run some parts of the model calculation in parallel. The number of cores that will be used can be defined via the \code{ncores} argument in the model input (see \code{\link{model_options}} or \code{\link{runbLS}}). Parallel computing is done using the package \code{snow} (and parts of \code{\link{snowfall}} and \code{\link{parallel}}). If the model is run in parallel mode, two parts of the model will be processed parallel, a) the core C++ routine, calculating the trajectories and corresponding touchdowns, and b) the calculation of the model results (i.e. \emph{C/E} etc.) from the touchdown catalogs.
}

\section{To Do:}{
\itemize{
\item{Plots (Contour, etc)}
\item{explain output csv}
\item{Docu}
\item{data.table}
}
}

\references{
Flesch, T. K., J. D. Wilson, et al. (2004). \dQuote{Deducing ground-to-air emissions from observed trace gas concentrations: A field trial.} Journal of Applied Meteorology 43(3): 487-502.

Flesch, T. K. (1996). \dQuote{The footprint for flux measurements, from backward Lagrangian stochastic models.} Boundary-Layer Meteorology 78(3-4): 399-404.

Flesch, T. K., J. D. Wilson, et al. (1995). \dQuote{Backward-time Lagrangian stochastic dispersion models and their application to estimate gaseous emissions.} Journal of Applied Meteorology 34(6): 1320-1332.
}

\keyword{ package }
\seealso{
	\code{\link{runbLS}}, \code{\link{genInputList}}, \code{\link{genSensors}}, \code{\link{genSources}}, \code{\link{genInterval}}, \code{\link{model_options}}, \code{\link{coreModel}}.
}
\examples{
\dontrun{
## set up some sensors and sources:
Sensors <- genSensors(
  PointSensor = list(x=0,y=0,z=2)
  ,LineSensor = list(x=c(-10,10),y=0,z=2,d=0.5)
  )
Source <- genSources(Source1=list('c',M=c(0,20),R=10))

## get default Interval data with optimized MaxFetch:
Int <- genInterval(MaxFetch=-20)

## generate TD catalog with temporary touchdown catalogs:
Model <- model_options(TDwrite=FALSE)
InList <- genInputList(Sensors,Source,Int,Model)
Run <- runbLS(InList,Cat.Path=getwd())

## look at catalog
TDs <- readCatalog(getCatalogs(Run,rn=1,Sensor='PointSensor')[,Catalog])	
plot(TDs,pch=20)

## temporary TD catalogs will be deleted on exiting R or
## by running the funtion cleanTemporary()
}
}
