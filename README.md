# bLSmodelR - An R package to set up and run a short-range atmospheric dispersion model

<!-- badges: start -->
[![R-CMD-check](https://github.com/ChHaeni/bLSmodelR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ChHaeni/bLSmodelR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Description
The package bLSmodelR provides functions to run a backward Lagrangian stochastic (bLS) dispersion model with the [R programming language](https://www.r-project.org/). Results from bLSmodelR are consistent with results from the freeware [WindTrax](http://www.thunderbeachscientific.com/).

## Installation

### Install package dependencies
```r
# install required packages
install.packages(c('Rcpp', 'rlecuyer', 'data.table', 'qs'))

# install optional packages for footprint plotting
# install.packages(c('sp', 'rgeos', 'geosphere', 'RgoogleMaps', 'maptools'))
```

### Install package from source
```r
# install bLSmodelR
devtools::install_github('ChHaeni/bLSmodelR', dependencies = FALSE)
```

## Bugs and Contact

If you have questions related to the package, it is best to contact me by email at contact@blsmodelr.slmail.me

Bugs and issues should be reported either through GitHub or via email at bugs@blsmodelr.slmail.me

## How to run the model
Guide to bLSmodelR: [https://github.com/ChHaeni/bLSmodelR/blob/main/Guide2bLSmodelR.r](https://github.com/ChHaeni/bLSmodelR/blob/main/Guide2bLSmodelR.r)

## bLS model
The model is a first-order Lagrangian stochastic model that is run in backward mode (i.e. backward in time) assuming horizontally homogeneous and vertically inhomogeneous Gaussian turbulence. The model is based on the paper published by [Flesch et al. (2004)](#Fl04). The basis of the model calculation is a generalized Langevin equation given as:

*du<sub>i</sub> = a<sub>i</sub>(x, u)\*dt + b<sub>ij</sub>(x, u)\*dξ<sub>j</sub>*
*dx<sub>i</sub> = u<sub>i</sub>\*dt*

where *dξ* is a random increment from a Gaussian distribution *N(0, dt)*. The indices *i, j = 1, 2, 3* represent the alongwind, crosswind and vertical axis of the rotated coordinate system. *x<sub>i</sub>* and *u<sub>i</sub>* are position and velocity of the trajectory at time *t*. An ensemble of upwind trajectories is calculated starting from the sensor position *x = (0, 0, z<sub>meas</sub> - d)* and the position and velocity of trajectory touchdowns on ground (i.e. on a horizontal plane at *x<sub>3</sub> = z<sub>0</sub> + d*, where the modelled, average wind speed equals 0 m/s.) are recorded to calculate the sensor concentration to source strength relationship as the sum of the inverse touchdown velocities, summed over all touchdowns within the source:

*C/E = 2/N<sub>traj</sub>\*Σ<sub>insideSource</sub>(1/w<sub>TD</sub>)*

where *C* represents the sensor concentration in ***mass m<sup>-3</sup>*** and source areas are defined as homogenously emitting sources, having emission rate *E* in ***mass m<sup>-2</sup> s<sup>-1</sup>***, since the *C/E* relationship is expressed in ***s m<sup>-1</sup>***. Analogously, the flux to source strength ratios (*<u'C'>/E*, *<v'C'>/E*, *<w'C'>/E*) are calculated, mainly for completeness and for possible vertical flux footprint corrections. Trajectories touching ground are perfectly reflected and their velocities are reversed to maintain covariances past reflection. In the main model, no deposition is modelled when trajectory touchdown happens. Furthermore, this model is not capable of modelling chemistry, and trace gases have to be assumed to be inert within the relevant travelling time. For further details on the model and the concentration to emission rate relationship, see [Flesch et al. (2004)](#Fl04), [Flesch (1996)](#Fl96) and partly [Flesch et al. (1995)](#Fl95).

## References
<a name="Fl95"></a>Flesch, T. K., J. D. Wilson, et al. (1995). "Backward-time Lagrangian stochastic dispersion models and their application to estimate gaseous emissions." Journal of Applied Meteorology 34(6): 1320-1332.

<a name="Fl96"></a>Flesch, T. K. (1996). "The footprint for flux measurements, from backward Lagrangian stochastic models." Boundary-Layer Meteorology 78(3-4): 399-404.

<a name="Fl04"></a>Flesch, T. K., J. D. Wilson, et al. (2004). "Deducing ground-to-air emissions from observed trace gas concentrations: A field trial." Journal of Applied Meteorology 43(3): 487-502.
