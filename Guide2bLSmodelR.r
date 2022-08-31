
#################################
# Guide to bLSmodelR Version 4+
#################################

# I am very sorry for the underdocumented help pages in bLSmodelR.
# I hope this guide will bring some better insight into the bLSmodelR...
# The current version of bLSmodelR to use with this guide is version 4.3-0 (2020-06-12)

library(bLSmodelR)

####################################
#~~~~~~~I. Setting up a Run:~~~~~~~#
####################################

#******
# I. What input is needed?
#******

	#******
	# I.M. What input is mandatory?

	#******
	# I.M.0: Have a look at the main model function *runbLS*:
		# ?runbLS
		# to run bLSmodelR a list of data.frames defining the model input is needed
		# ?genInputList
		# The argument 'ModelInput' in *runbLS* contains 2 mandatory entries:
		# $Interval 	-> I.M.1
		# $Sensors 		-> I.M.2
		#
		# and 3 optional entry:
		# $Sources 		-> I.O.1
		# $Model 		-> I.O.2
		# $Tolerances 	-> I.O.3
	#******

	#******
	# I.M.1: $Interval:
		# ?genInterval
		print(Interval <- genInterval())
		str(Interval)
		# contains all information about the individual intervals that will be calculated
		# i.e. necessary turbulence and windprofile parameters incl. wind direction
		# Sensors and Sources that will be included in the calculation of the interval can be defined as well
		# if arguments 'SensorNames' and/or 'SourceNames' are not provided, all available Sensors/Sources will be included
		# when gathering the model input with the function *genInputList* (c.f. example below IntervalDemo).
		# example 1: 'Sensor1' will be calculated in first interval, 'Sensor2' in second
		print(Interval1 <- genInterval(SensorNames = c("Sensor1","Sensor2")))
		# example 2: 'Sensor1' and 'Sensor2' will be calculated in the same interval
		print(Interval2 <- genInterval(SensorNames = "Sensor1,Sensor2"))
		### 'MaxFetch' argument
		# argument 'MaxFetch' can be provided as negative number. If done so, the maximum distance between the specific 
		# Sources and the Sensors will be calculated based on the wind direction during the model run
		# and the absolute value of the 'MaxFetch' argument will be added on top of this distance in order to define the MaxFetch 
		# values (i.e. the maximum distance upwind of the Sensor) for each interval (see example in ?genInterval)
		print(Interval3 <- genInterval(MaxFetch = -50))
		####
		# an option to construct the Interval data.frame by hand is given as well (using argument 'Data').
		# This comes most of the time handier and additionally allows for more information changing on an interval basis (e.g. start/end times etc.), 
		# which is passed to the output. Columns that should be recognized as Interval arguments should be named according to the 
		# corresponding arguments, i.e names(as.list(args(genInterval)))[2:14]
		IntervalData <- cbind(Interval3,Numeric=c(1.2,3.1),Characters=letters[1:2],Factors=factor(letters[1:2]),POSXct=as.POSIXct(1:2,origin=Sys.time()),stringsAsFactors=FALSE)
		print(IntervalDemo <- genInterval(Data=IntervalData))
		str(IntervalDemo)
	#******

	#******
	# I.M.2: $Sensors:
		# ?genSensors
		# genSensors accepts either lists or data.frames as argument input
		# these lists/data.frames need specific named entries (mandatory are
		# entries x, y, and z; optional is entry name - the name of the sensor):
		# Point Sensors need 3 entries (x,y,z)
		# PathSensors need one of these entries (x,y,z) with length >1 and can 
		# have two additional arguments:
		# 	- entry 'd' refines the way to interpolate the path sensor by multiple point sensors:
		# 		-> d = Distance between point sensors in m
		#		- entry 'id' defines 'sub-groups' of a sensor, such that multiple sensors at different 
		#			locations can be treated as one sensor 
		# providing no argument to genSensors results in an error
		# genSensors() 
		print(Sensors <- genSensors(PointSensor=list(x=0,y=0,z=1.2)))
		str(Sensors)
		# other example:
		print(Sensors1 <- genSensors(data.frame(name=paste0("Sensor",1:3),x=10,y=0,z=c(0.5,1.2,2.05))))
		# defining Path Sensors:
		print(Sensors2 <- genSensors(
			LineSensor1=list(x=c(0,20),y=0,z=1.2),
			LineSensor2=list(x=c(-10,20),y=10,z=c(2,2.2),d=0.5),
			PathSensor1=list(x=c(0, 10, 10),y=c(15,25,30),z=1.2),
			PathSensor2=list(x=c(2, 2, 10, 10, 20),y=c(0, 15,5,10, 5),z=1.2,id=c(1,1,2,2,2))
		))
		plot(Sensors2[Sensors2$"Sensor Name" == "PathSensor2",])
		# joining Sensors:
		head(join(Sensors1,Sensors2))

		# Now define the two Sensors included in IntervalDemo for further use:
		print(SensorsDemo <- genSensors(Sensor1=list(x=10,y=0,z=1.2),Sensor2=list(x=20,y=0,z=2.05)))
	#******

	#****** 
	# Gathering ModelInput:
		# ?genInputList
		# Given these 2 list entries, bLSmodelR could be run to pre-calculate TD catalogs - but only, if TDonly==TRUE:
		# *genInputList* will set TDonly==TRUE with a warning
		# Default values for Model Parameters and Tolerances (see I.O.2 and I.O.3 below) will be used
		InputWithoutSource <- genInputList(IntervalDemo,SensorsDemo)
		print(InputWithoutSource)
		# Since no sensors were defined explicitly in the Interval data.frame, all existing sensors were added to each interval.
	#****** 




	#******
	# I.O. Optional model input:
	#******
	# I.O.1: $Sources:
		# Sources are treated as homogeneously emitting source ares, defined by polygons.
		# ?genSources
		print(Sources <- genSources(list(x=rnorm(10),y=rnorm(10))))
		str(Sources)
		# generating sources is done in a somewhat inconvenient way...
		# source areas are either defined by: 
		# 1) a list of named lists, where the names of the list entries correspond to the source area names.
		# 2) arguments providing lists or data.frames, where the argument names correspond to the source are names.
		# There are two optional polygon "Classes" to help building source areas:
		# circle 'Class'="c": 
		print(Circle <- genSources(SourceCirc=list(Class="c",M=c(10,50),R=10)))
		# rectangle 'Class'="r": 
		print(Rect <- genSources(SourceRect=list(Class="r",x1=0,x2=10,y1=0,y2=20)))
		# Otherwise source polygons are defined by a list with vectors x and y: 
		print(Poly <- genSources(SourcePoly=list(
			x=c(seq(-10,10,1),seq(9.5,-9.5,-0.5))-20,
			y=c(rep(0,21),seq(h <- sqrt(20^2-10^2)/20,h*20,h),seq(19*h,h,-h))+30
			)))
		# or equivalently by supplying a data.frame:
		print(Poly <- genSources(SourcePoly=data.frame(
			x=c(seq(-10,10,1),seq(9.5,-9.5,-0.5))-20,
			y=c(rep(0,21),seq(h <- sqrt(20^2-10^2)/20,h*20,h),seq(19*h,h,-h))+30
			)))
		# Multiple entries with equal names are treated as one source.
		# One can provide them as one list argument with named list entries
		print(MultiPoly <- genSources(list(
			SourceRects=list(Class="r",x1=10,x2=20,y1=5,y2=25),
			SourceRects=list(Class="r",x1=-10,x2=-20,y1=5,y2=25),
			SourceRects=list(Class="r",x1=-10,x2=10,y1=-10,y2=0),
			OtherSource=list(Class="r",x1=-5,x2=5,y1=20,y2=30)			
			)))
		# or, equivalently, as multiple (possibly named) arguments:
		print(MultiPoly <- genSources(
			SourceRects=list(Class="r",x1=10,x2=20,y1=5,y2=25),
			SourceRects=list(Class="r",x1=-10,x2=-20,y1=5,y2=25),
			SourceRects=list(Class="r",x1=-10,x2=10,y1=-10,y2=0),
			OtherSource=list(Class="r",x1=-5,x2=5,y1=20,y2=30)			
			))
		# join all previous in one for demo:
		print(SourcesDemo <- join(
			Circle,
			Rect,
			Poly,
			MultiPoly
			))
		head(SourcesDemo)

	#******


	#******
	# I.O.2: $Model:
		# ?genModel
		print(Model <- genModel())
		str(Model)
		# kv:			von Karman constant
		# A:			Scaling Factor (see Flesch et al., 2004)
		# alpha:		fraction of "tau_L" to define timestep (see Flesch et al., 2004)
		# wTDcutoff:	lower limit of touchdown velocities (TD velocities below the cutoff will be set to the cutoff value)		
		# TDwrite:		Writing (i.e. permanently saving) touchdown catalogs?		
		# overwriteTD:	Overwrite existing (identical) touchdown catalogs?
		# TDread:		Reading previously saved touchdown catalogs?
		# TDonly:		Calculate touchdown catalogs only?
		# ncores:		Number of cores to use. 'ncores' <= 0 == "use all detected cores" as given by the 
		#				default call of *parallel::detectCores* (although this is not recommended 
		#				-> see: ?parallel::detectCores)
		#
		#  For demo: do not save TD catalogs:
		print(ModelDemo <- genModel(TDwrite=FALSE))
	#******

	#******
	# I.O.3: $Tolerances:
		# ?genTolerances
		print(Tol <- genTolerances())
		str(Tol)
		# only used for choosing appropriate touchdown catalogs, if TDread is "TRUE".
		# Tol.Zero prohibits using tolerances when choosing TD cataolgs,
		# i.e. only exactly matching catalogs will be used:
		print(genTolerances(Tol.Zero=TRUE))
	#******

	#****** 
	# Gathering ModelInput:
		# Generate ModelInput for demo:
		print(demoInput <- genInputList(IntervalDemo,SensorsDemo,SourcesDemo,ModelDemo,Tol))
		# If you use default values (as here Tolerances), you don't need to supply them:
		identicalDemoInput <- genInputList(IntervalDemo,SensorsDemo,SourcesDemo,ModelDemo)
		identical(demoInput,identicalDemoInput)
	#****** 


#******
# II. Have a look at your set up
#******

	#****** 
		siteMap(demoInput,polygon.args = list(col = rainbow(5,s=0.3)),
			sources.text.args = list(col = rainbow(5,v=0.5)))
		# or equivalently: 
		plot(demoInput,polygon.args = list(col = rainbow(5,s=0.3)),
			sources.text.args = list(col = rainbow(5,v=0.5)))
		args(siteMap)
		# additional possible calls (arguments order doesn't matter):
		par(mfrow=c(2,2))
		plot(SourcesDemo,SensorsDemo,main="plot(SourcesDemo,SensorsDemo)")
		plot(SensorsDemo,SourcesDemo,main="plot(SensorsDemo,SourcesDemo)")
		plot(SensorsDemo,main="plot(SensorsDemo)")
		plot(SourcesDemo,main="plot(SourcesDemo)")
		# make a somewhat nicer looking plot :-)
		nicerMap <- function(...){}
		body(nicerMap) <- expression({plot(demoInput
			,polygon.args=list(col=rainbow(5,s=0.3))
			,points.args=list(pch=4,cex=1,col=c("darkgreen","darkred"),lwd=2)
			,sensors.text.args=list(pos=4)
			,sources.text.args=list(col = rainbow(5,v=0.5))
			,...)})
		nicerMap()
		# add a scale bar:
		bnicerMap <- body(nicerMap)
		bnicerMap[3] <- expression(addScaleBar())
		body(nicerMap) <- bnicerMap
		nicerMap()
		bnicerMap[3] <- expression(addScaleBar(pos=4))
		body(nicerMap) <- bnicerMap
		nicerMap()
		# wind direction from NW would be best:
		bnicerMap[4] <- expression(addWindrose(WD=315,pos=1))
		body(nicerMap) <- bnicerMap
		nicerMap()
	#****** 


#******
# III. Run the model
#******

	#****** 
	# ?runbLS
		# define a temporary catalog path. Even if you don't save the TD catalogs at the end,
		# temporary TD catalogs will be stored there and later removed.
		# Make sure to have enough space available there. The following example will only require
		# around 2 MB for the temporary TD catalogs.
		Cat.Path <- getwd()
		
		# set Interval WD to 315 deg N for demo:
		demoInput$Interval[,"WD [deg N]"] <- 315

		# run the model
		args(runbLS)
		# results as data.table:
		DemoOutput_DataTable <- runbLS(demoInput,Cat.Path=Cat.Path)
		# temporary TD catalogs will be used as long as they're not removed (i.e. actively removed
		# by a call to cleanTemporary(), or by exiting the R session (where cleanTemporary() will be called before exiting))
		# results as data.frame
		DemoOutput <- runbLS(demoInput,Cat.Path=Cat.Path,asDT=FALSE)
		#
		DemoOutput
		head(DemoOutput,1)
		names(attributes(DemoOutput))
		# initially attached columns:
		DemoOutput[,39:42]
		str(DemoOutput[,39:42])
		# change from data.frame to data.table:
		setDT(DemoOutput)
		DemoOutput
		# and back to data.frame:
		setDF(DemoOutput)
		DemoOutput
		# switch column names to simple:
		DemoOutput <- switchNames(DemoOutput)
		DemoOutput
		# change back to data.table with more 'informative' column names:
		switchNames(setDT(DemoOutput),FALSE)
		DemoOutput
		# switch names back:
		switchNames(DemoOutput)

		# switch the heights of sensors (in a somewhat unconventional way). One catalog needs to be calculated
		# due to the new, longer Fetch required. This could have been omitted by setting a large enough MaxFetch value. 
		demoInput2 <- demoInput
		demoInput2$Sensors <- genSensors(data.frame(name=c("S1.High","S2.Low"),x=(1:2)*10,y=0,z=c(2.05,1.2)))
		demoInput2$Interval[,"Sensor Names (sep = \",\")"] <- c("S1.High","S2.Low")
		DemoOutput2 <- runbLS(demoInput2,Cat.Path=Cat.Path)
		DemoOutput2[,.(Sensor,Source,SensorHeight,CE)]
		DemoOutput[,.(Sensor,Source,SensorHeight,CE)]

		# extract Sensor1 and sort C/E in decreasing order (have a look at data.table::setorder for increasing/decreasing ordering -> i.e. -/+ for decreasing/increasing column ordering):
		S1Output <- extractResult(DemoOutput,sortArgs = list("Sensor"="Sensor1",-"CE"),dropAttr=FALSE)
		S1Output[]
		names(attributes(S1Output))
		# otherwise, if dropAttr=TRUE (default), model specific attributes are dropped:
		S1Results <- extractResult(DemoOutput,sortArgs = list("Sensor"="Sensor1",-"CE"))
		S1Results[]
		names(attributes(S1Results))
		#
		# extract Sources c("SourceRect","SourceRects") and drop all columns but c("Sensor","Source","C/E","C/E SE"):
		OutShort <- extractResult(DemoOutput,sortArgs = list("Source"=c("SourceRect","SourceRects")),keep=c("Sensor","Source","CE","CE_se"))
		OutShort
		#
		# get catalog: (somewhat inconvenient... will be improved someday)
		# row 1:
		CatInfo <- getCatalogs(DemoOutput,1)
		Catalog <- readCatalog(CatInfo[,Catalog])
		Catalog
		par(mfrow=c(1,1))
		plot(Catalog, col = "#00000011", pch = 20)
		
		# Visualize footprint (quite a number of arguments, sorry for that!):
		args(plotFootprint)
		par(mfrow=c(2,2))
		# Plot C-Footprint of Sensor 1:
		nicerMap(main="C-FP as modelled")
		plotFootprint(DemoOutput,"Sensor1",1,add=TRUE,addSource=FALSE,showSensor=FALSE,lpos="topright",showPerc=TRUE,alpha=0.5)
		# Plot 'despiked' C-Footprint of Sensor 1:
		nicerMap(main="C-FP (use.avg = TRUE)")
		plotFootprint(DemoOutput,"Sensor1",1,use.avg=TRUE,add=TRUE,addSource=FALSE,showSensor=FALSE,lpos="topright",showPerc=TRUE,alpha=0.5)
		# Plot 'despiked' C-Footprint of Sensor 1 using symmetry property due to horizontal homogeneity of turbulence:
		nicerMap(main="C-FP (use.avg = TRUE, use.sym = TRUE)")
		plotFootprint(DemoOutput,"Sensor1",1,use.avg=TRUE,use.sym=TRUE,add=TRUE,addSource=FALSE,showSensor=FALSE,lpos="topright",showPerc=TRUE,alpha=0.5)
		
		# Plot w'C'-Footprint:
		x11(width=10,height=5)
		par(mfrow=c(1,2))
		nicerMap(main="w'C'-FP")
		plotFootprint(DemoOutput,"Sensor1",1,type="wCE",use.avg=TRUE,add=TRUE,addSource=FALSE,showSensor=FALSE,lpos="topright",showPerc=TRUE)
		# Plot C-Footprint of Sensor 1 and Sensor 2:
		nicerMap(main="C-FP of both Sensors (units: CE s/m)")
		plotFootprint(DemoOutput,"Sensor1",1,add=TRUE,use.avg=TRUE,use.sym=TRUE,addSource=FALSE,showSensor=FALSE,lpos="topright",alpha=0.5,breaks = function(x) exp(seq(log(1E-4),log(0.027),length.out=4)),sigNums=2)
		plotFootprint(DemoOutput,"Sensor2",2,add=TRUE,use.avg=TRUE,use.sym=TRUE,addSource=FALSE,showSensor=FALSE,lpos=NA,alpha=0.3,breaks = function(x) exp(seq(log(1E-4),log(0.027),length.out=4)))

	#****** 

#******
# IV. Deposition post-processing
#******
		# Run the deposition post-processing function: deposition
		args(deposition)
		# vDep must be given as m/s. It is defined as the 'surface' deposition velocity ath height d + z0 (the "bLS model surface").
		# post-process results with a 'rather high' deposition velocity of 3 cm/s: 
		RunDep <- deposition(DemoOutput,0.03,Sensor="Sensor1",Source="SourceRects")
		RunDep
		RunDep[,.(CE,CE_Dep,"Reduction of CE by deposition (in %)"=(CE - CE_Dep)/CE*100)]

		# spatially inhomogeneous deposition velocity 
		# this is still experimental and input format may change in future updates...
		vDepAreas <- join(Circle,Rect,MultiPoly)
		vDepList <- list(
			SourceRects = 0
			,SourceRect = 0.03
			,OtherSource = 0.01
			)
		plot(join(Poly,vDepAreas),SensorsDemo,polygon.args = list(col = rainbow(5,s=0.3)),sources.text.args = list(
			labels = paste0("vDep = ",c("0 (emitting area)","0.002 (same as default,\n because it is not defined in \"vDepList\")","0.03","0","0.01"))
			,col = rainbow(5,v=0.5)))
		text(0,35,"vDep = 0.002 (all non-defined areas)",cex=0.5,font=2)
		RunDep2 <- deposition(DemoOutput,0.002,Sensor="Sensor1",Source="SourcePoly",vDepSpatial = list(vDepList,vDepAreas))
		RunDep2[,.(CE,CE_Dep,"Reduction of CE by deposition (in %)"=(CE - CE_Dep)/CE*100)]

        #### example on new vDepSpatial treatment
        ## define key
        #snames <- unique(Sources[, 1])
        ## source names key
        #vdn_key <- setNames(snames, snames[c(2, 1, 4, 3)])
        ## 

        ## 1. data.frame-like (is.data.frame/inherits, one column defining patches/polygons/"Sources" per row/interval, columns with patches names defining vDep per row/interval
        #vds <- nodep[, .(rn, Sensor, Source)]
        #vds[, Spatial := nodep[, vdn_key[Source]]]
        ## spatial, zones, sources, select, regions
        #vds[, id := .I]
        #vdadd <- dcast(vds, id + rn + Sensor + Source + Spatial ~ Spatial, value.var = 'Spatial', fun.aggregate = function(x) 0.0)
        ## check:
        ##   - empty row entry
        #vdadd[sample.int(.N, 10), Spatial := '']

        ## define vDep via columns & set vDep to either vDep or 0
        #vDepSpatialList <- list(
        #    # list with vDep == 0
        #    vdadd,
        #    # Sources object
        #    Sources
        #) 

        ## run deposition
        #dep <- deposition(nodep, vDep = 'vDep', vDepSpatial = vDepSpatialList, ncores = 1)


#******
# V. Parallelism
#******
		# A very! simple attempt to parallelize parts of the model calculation has been implemented using the package snow (and snowfall)
		# It is possible to either: 
		# a) provide the number of cores (argument: ncores) to build a SOCKET type structure or 
		DemoOutput_Parallel <- runbLS(demoInput,Cat.Path=Cat.Path,ncores=2)
		# b) to have a network already running before calling runbLS (if so, keep in mind, that the TD catalog must be available to all nodes!!!)
		sfInit(TRUE,2)
		DemoOutput_Parallel2 <- runbLS(demoInput,Cat.Path=Cat.Path);sfStop()

		
