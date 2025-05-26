\name{Undocumented Functions}
\alias{Undocumented}
\alias{\%w/o\%}
\alias{addScaleBar}
\alias{addWindrose}
\alias{assignCat}
\alias{check.compact}
\alias{cleanTemporary}
\alias{compactCatalog}
\alias{extractResult}
\alias{getArea}
\alias{getCatalogs}
\alias{initializeCatalog}
\alias{join}
\alias{join.bLSresult}
\alias{join.deposition}
\alias{prepareIntervals}
\alias{plot.InputList}
\alias{plotFootprint}
\alias{readCatalog}
\alias{rotateCatalog}
\alias{rotate}
\alias{round125}
\alias{siteMap}
\alias{setDT}
\alias{setDF}
\alias{sortData}
\alias{tag_inside}
\alias{uvw0}
\alias{updatePath}
\alias{writeCatalog}
\alias{rebuildCatListFile}
\alias{switchNames}
\alias{CalcSteps}
\alias{CatPath}
\alias{Catalogs}
\alias{ModelInput}
\alias{ConcPalette}
\alias{FluxPalette}
\alias{add_data}
\alias{avgCE_sensors}
\alias{avgCE_sensors_se}
\alias{avgCE_sources}
\alias{avgCE_sources_se}
\alias{avgSourcesWeights}
\alias{calc_fct}
\alias{combineSensors}
\alias{combineSources}
\alias{convert}
\alias{getPathLengths}
\alias{inters}
\alias{plot.bLSresult}
\alias{plotPolys}
\alias{pnt.in.poly}
\alias{print.InputList}
\alias{print.Interval}
\alias{splitSource}
\alias{collect_results}
\alias{depoSlurm}
\alias{find_partition}
\alias{rename_sensor}
\alias{rename_source}
\alias{runbLSlurm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
undocumented functions
}
\description{
List of undocumented bLSmodelR functions:
}

\details{
\itemize{
\item{\code{\%w/o\%}: }{x without y, analogous to \code{\link{\%in\%}}}
\item{\code{addScaleBar}: }{add a scale bar to the current plot}
\item{\code{addWindrose}: }{add a wind rose (arrow showing wind direction) to the current plot}
\item{\code{assignCat}: }{assign \code{\link{coreModel}} output to an empty, initialized TD catalog}
\item{\code{CalcSteps}: }{Get attribute CalcSteps from bLS output}
\item{\code{Catalogs}: }{Get attribute Catalogs from bLS output}
\item{\code{CatPath}: }{Get attribute CatPath from bLS output}
\item{\code{check.compact}: }{check TD catalog for compact attributes}
\item{\code{cleanTemporary}: }{remove temporary TD catalogs}
\item{\code{compactCatalog}: }{make TD catalog compact for saving (not reversible due to rounding!)}
\item{\code{ConcPalette}: }{Colour palette}
\item{\code{extractResult}: }{extract results from model output}
\item{\code{FluxPalette}: }{Colour palette}
\item{\code{getArea}: }{get area of a polygon}
\item{\code{getCatalogs}: }{get TD catalog pathnames from model output}
\item{\code{initializeCatalog}: }{initialize TD catalog}
\item{\code{join}: }{see \code{\link{join.Sensors}} and \code{\link{join.Sources}}}
\item{\code{join.bLSresult}: }{merge bLS results consistently}
\item{\code{join.deposition}: }{merge deposition results consistently}
\item{\code{ModelInput}: }{Get attribute ModelInput from bLS output}
\item{\code{prepareIntervals}: }{link existing and upcoming TD catalogs to intervals and possibly optimize MaxFetch}
\item{\code{plot.InputList}: }{S3method calling \code{siteMap}}
\item{\code{plotFootprint}: }{plot C-/u'C'-/w'C'-Footprint}
\item{\code{readCatalog}: }{read TD catalog}
\item{\code{rebuildCatListFile}: }{rebuild file with TD catalog list}
\item{\code{rotate}: }{rotate coordinates around a fix point}
\item{\code{rotateCatalog}: }{rotate Catalog (into WD)}
\item{\code{round125}: }{round digits to nearest predefined numbers}
\item{\code{siteMap}: }{plot site map}
\item{\code{sortData}: }{sort data, partially adopted from \code{\link[data.table:setorderv]{setorderv}}}
\item{\code{tag_inside}: }{tag TDs inside source area}
\item{\code{updatePath}: }{update path to TD catalogs}
\item{\code{uvw0}: }{get initial (t=0) velocities from TD catalog}
\item{\code{writeCatalog}: }{write TD catalog (default: saved as compact TD catalog)}
\item{\code{switchNames}: }{change column names (simple/descriptive) of model output}
\item{\code{setDT}: }{change model output to a \code{data.table}, see \code{\link[data.table:setDT]{setDT}}}
\item{\code{setDF}: }{change model output to a \code{data.frame}, see \code{\link[data.table:setDF]{setDF}}}
\item{\code{add_data}}
\item{\code{avgCE_sensors}}
\item{\code{avgCE_sensors_se}}
\item{\code{avgCE_sources}}
\item{\code{avgCE_sources_se}}
\item{\code{avgSourcesWeights}}
\item{\code{calc_fct}}
\item{\code{combineSensors}}
\item{\code{combineSources}}
\item{\code{convert}}
\item{\code{getPathLengths}}
\item{\code{inters}}
\item{\code{plot.bLSresult}}
\item{\code{plotPolys}}
\item{\code{pnt.in.poly}}
\item{\code{print.InputList}}
\item{\code{print.Interval}}
\item{\code{splitSource}}
\item{\code{collect_results}}
\item{\code{depoSlurm}}
\item{\code{find_partition}}
\item{\code{rename_sensor}}
\item{\code{rename_source}}
\item{\code{runbLSlurm}}

}
}


\author{
Christoph Haeni
}
\note{
some of these functions will be documented or dropped in future package versions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bLSmodelR-package}}, \code{\link{runbLS}}
}
