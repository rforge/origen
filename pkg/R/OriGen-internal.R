#This is R code for OriGen made by John Michael O. Ranola ranolaj@uw.edu
#if the function is not to be accessed by users, start with a period(.)

.is.wholenumber <-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol



ConvertPEDData<-function(PlinkFileName,LocationFileName){
#DataFileName should be the base name of plink ped/map format
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#Location file should be space/tab delimited with columns ID,Label,AltLabel,Long,Lat
	NumberSNPs=1
	MapFileName=paste(PlinkFileName,".map",sep="")
	temp=.Fortran("COUNT_NUMBER_LINES",NumberSNPs=as.integer(NumberSNPs),PlinkFileName=as.character(MapFileName),PACKAGE="OriGen")
	NumberSNPs=temp$NumberSNPs
	print(c("NumberSNPs",NumberSNPs))
	print(MapFileName)
	
	SampleSites=1
	temp=.Fortran("COUNT_NUMBER_SAMPLE_SITES",SampleSites=as.integer(SampleSites),LocationFileName=as.character(LocationFileName),PACKAGE="OriGen")
	SampleSites=temp$SampleSites
	print(SampleSites)
	
	DataArray=array((1:(2*SampleSites*NumberSNPs)),c(2,SampleSites,NumberSNPs))
	DataArray=array(0,c(2,SampleSites,NumberSNPs))
	SampleCoordinates=array(0,c(SampleSites,2))
	MembersList=array('xx',SampleSites)
	
#the following includes member names in the code, however r cannot pass vectors of strings to fortran... need a fix...	#ResultsRaw=.Fortran("FORMAT_PLINK_DATA",DataArray=as.integer(DataArray),SampleCoordinates=as.double(SampleCoordinates),MembersList=as.character(MembersList),PlinkFileName=as.character(paste(PlinkFileName,".ped",sep="")),LocationFileName=as.character(LocationFileName),SampleSites=as.integer(SampleSites),NumberSNPs=as.integer(NumberSNPs))
	
ResultsRaw=.Fortran("FORMAT_PLINK_DATA",DataArray=as.integer(DataArray),SampleCoordinates=as.double(SampleCoordinates),PlinkFileName=as.character(paste(PlinkFileName,".ped",sep="")),LocationFileName=as.character(LocationFileName),SampleSites=as.integer(SampleSites),NumberSNPs=as.integer(NumberSNPs),PACKAGE="OriGen")
	#ResultsRaw=.Fortran("RTRIALS",DataArray=as.integer(DataArray),DataLength=as.integer(c(2,SampleSites,NumberSNPs)))
	
ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
	
	return(ResultsRaw)
}


FitOriGenModel<-function(DataArray,SampleCoordinates,MaxGridLength=20,RhoParameter=10){
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#SampleCoordinates[SampleSites,2] gives the locations of the grouped data
#This function takes in the data, fits the model, and returns the allele frequency surfaces
	if(!.is.wholenumber(MaxGridLength)){
		stop("MaxGridLength must be an integer")
	}
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(RhoParameter<=0){
		stop("RhoParameter must be greater than 0")
	}
	if(length(SampleCoordinates[1,])!=2){
		stop("SampleCoordinates should give the Long/Lat coordinates of the grouped data so it should only contain 2 columns")
	}
	NumberSNPs=length(DataArray[1,1,])
	SampleSites=length(DataArray[1,,1])
	GridLength=array(0,2)
	GridCoordinates=array(0.,dim=c(2,MaxGridLength))
	
	GridAndCoordResults=.Fortran("UPDATE_GRID_COORD_SQUARE2",GridCoordinates=as.double(GridCoordinates),SampleCoordinates=as.double(SampleCoordinates),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),SampleSites=as.integer(SampleSites),PACKAGE="OriGen")
	
	GridLength=GridAndCoordResults$GridLength
	GridCoordinates=GridAndCoordResults$GridCoordinates
	
	AlleleFrequencySurfaces=array(0,dim=c(NumberSNPs,GridLength[1],GridLength[2]))
	ResultsRaw=.Fortran("FITORIGENMODEL",AlleleFrequencySurfaces=as.double(AlleleFrequencySurfaces),DataArray=as.integer(DataArray),NumberSNPs=as.integer(NumberSNPs),GridLength=as.integer(GridLength),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),MaxGridLength=as.integer(MaxGridLength),SampleCoordinates=as.double(SampleCoordinates),GridCoordinates=as.double(GridCoordinates),PACKAGE="OriGen")
	
	ResultsRaw$AlleleFrequencySurfaces=array(ResultsRaw$AlleleFrequencySurfaces,c(NumberSNPs,GridLength[1],GridLength[2]))
	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
ResultsRaw$GridCoordinates=array(ResultsRaw$GridCoordinates,c(2,MaxGridLength))
	
	return(ResultsRaw)
}





ConvertUnknownPEDData<-function(PlinkFileName,LocationFileName,PlinkUnknownFileName){
#DataFileName should be the base name of plink ped/map format
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#Location file should be space/tab delimited with columns ID,Label,AltLabel,Long,Lat

	NumberSNPs=1
temp=.Fortran("COUNT_NUMBER_LINES",NumberSNPs=as.integer(NumberSNPs),PlinkFileName=as.character(paste(PlinkFileName,".map",sep="")),PACKAGE="OriGen")
	NumberSNPs=temp$NumberSNPs
	print(c("NumberSNPs",NumberSNPs))
	
	NumberKnown=1
temp=.Fortran("COUNT_NUMBER_LINES",NumberKnown=as.integer(NumberKnown),PlinkFileName=as.character(paste(PlinkFileName,".ped",sep="")),PACKAGE="OriGen")
	NumberKnown=temp$NumberKnown

	NumberUnknowns=1
temp=.Fortran("COUNT_NUMBER_LINES",NumberUnknowns=as.integer(NumberUnknowns),PlinkUnknownFileName=as.character(paste(PlinkUnknownFileName,".ped",sep="")),PACKAGE="OriGen")
	NumberUnknowns=temp$NumberUnknowns
	
	SampleSites=1
temp=.Fortran("COUNT_NUMBER_SAMPLE_SITES",SampleSites=as.integer(SampleSites),LocationFileName=as.character(LocationFileName),PACKAGE="OriGen")
	SampleSites=temp$SampleSites
	print(SampleSites)
	
	DataArray=array((1:(2*SampleSites*NumberSNPs)),c(2,SampleSites,NumberSNPs))
	DataArray=array(0,c(2,SampleSites,NumberSNPs))
	SampleCoordinates=array(0,c(SampleSites,2))
	MembersList=array('xx',SampleSites)
	Membership=1:NumberKnown
	
#the following includes member names in the code, however r cannot pass vectors of strings to fortran... need a fix...	
#ResultsRaw=.Fortran("FORMAT_PLINK_DATA",DataArray=as.integer(DataArray),SampleCoordinates=as.double(SampleCoordinates),MembersList=as.character(MembersList),PlinkFileName=as.character(paste(PlinkFileName,".ped",sep="")),LocationFileName=as.character(LocationFileName),SampleSites=as.integer(SampleSites),NumberSNPs=as.integer(NumberSNPs))
	
UnknownData=array(0,c(NumberUnknowns,NumberSNPs))
ResultsRaw=.Fortran("FORMAT_UNKNOWN_PLINK_DATA",DataArray=as.integer(DataArray),SampleCoordinates=as.double(SampleCoordinates),PlinkFileName=as.character(paste(PlinkFileName,".ped",sep="")),LocationFileName=as.character(LocationFileName),SampleSites=as.integer(SampleSites),NumberSNPs=as.integer(NumberSNPs),PlinkUnknownFileName=as.character(paste(PlinkUnknownFileName,".ped",sep="")),NumberUnknowns=as.integer(NumberUnknowns),UnknownData=as.integer(UnknownData),Membership=as.integer(Membership),NumberKnown=as.integer(NumberKnown),PACKAGE="OriGen")

ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
ResultsRaw$UnknownData=array(ResultsRaw$UnknownData,c(NumberUnknowns,NumberSNPs))
	
	return(ResultsRaw)
}



FitOriGenModelFindUnknowns<-function(DataArray,SampleCoordinates,UnknownData,MaxGridLength=20,RhoParameter=10){
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#by jmor
#SampleCoordinates[SampleSites,2] gives the locations of the grouped data
#UnknownData[NumberUnknowns,NumberSNPs] gives the number of major alleles for the current unknown individual
#This function takes in the data, fits the model, and returns the allele frequency surfaces
	if(!.is.wholenumber(MaxGridLength)){
		stop("MaxGridLength must be an integer")
	}
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(RhoParameter<=0){
		stop("RhoParameter must be greater than 0")
	}
	if(length(SampleCoordinates[1,])!=2){
		stop("SampleCoordinates should give the Long/Lat coordinates of the grouped data so it should only contain 2 columns")
	}
	NumberSNPs=length(DataArray[1,1,])
	SampleSites=length(DataArray[1,,1])
	GridLength=array(0,2)
	GridCoordinates=array(0.,dim=c(2,MaxGridLength))
	
	GridAndCoordResults=.Fortran("UPDATE_GRID_COORD_SQUARE2",GridCoordinates=as.double(GridCoordinates),SampleCoordinates=as.double(SampleCoordinates),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),SampleSites=as.integer(SampleSites),PACKAGE="OriGen")
	
	GridLength=GridAndCoordResults$GridLength
	GridCoordinates=GridAndCoordResults$GridCoordinates
	
	NumberUnknowns=length(UnknownData[,1])
	UnknownGrids=array(0,dim=c(GridLength[1],GridLength[2],NumberUnknowns))
	ResultsRaw=.Fortran("FITORIGENMODELFINDUNKNOWNS",UnknownGrids=as.double(UnknownGrids),DataArray=as.integer(DataArray),NumberSNPs=as.integer(NumberSNPs),GridLength=as.integer(GridLength),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),MaxGridLength=as.integer(MaxGridLength),SampleCoordinates=as.double(SampleCoordinates),NumberUnknowns=as.integer(NumberUnknowns),UnknownData=as.integer(UnknownData),GridCoordinates=as.double(GridCoordinates),PACKAGE="OriGen")
	
	ResultsRaw$UnknownGrids=array(ResultsRaw$UnknownGrids,c(GridLength[1],GridLength[2],NumberUnknowns))
	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
ResultsRaw$GridCoordinates=array(ResultsRaw$GridCoordinates,c(2,MaxGridLength))

ResultsRaw$UnknownData=array(ResultsRaw$UnknownData,c(NumberUnknowns,NumberSNPs))
	
	return(ResultsRaw)
}



FindRhoParameterCrossValidation<-function(PlinkFileName,LocationFileName,MaxIts=6,MaxGridLength=20){
#DataFileName should be the base name of plink ped/map format
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#Location file should be space/tab delimited with columns ID,Label,AltLabel,Long,Lat
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(MaxIts<3){
		stop("MaxIts must be greater than 3")
	}
	
	NumberSNPs=1
	MapFileName=paste(PlinkFileName,".map",sep="")
	PedFileName=paste(PlinkFileName,".ped",sep="")
	temp=.Fortran("COUNT_NUMBER_LINES",NumberSNPs=as.integer(NumberSNPs),PlinkFileName=as.character(MapFileName),PACKAGE="OriGen")
	NumberSNPs=temp$NumberSNPs
	print(c("NumberSNPs",NumberSNPs))
	print(MapFileName)
	
	SampleSites=1
	temp=.Fortran("COUNT_NUMBER_SAMPLE_SITES",SampleSites=as.integer(SampleSites),LocationFileName=as.character(LocationFileName),PACKAGE="OriGen")
	SampleSites=temp$SampleSites
	print(SampleSites)
	
	DataArray=array((1:(2*SampleSites*NumberSNPs)),c(2,SampleSites,NumberSNPs))
	DataArray=array(0,c(2,SampleSites,NumberSNPs))
	SampleCoordinates=array(0,c(SampleSites,2))
	MembersList=array('xx',SampleSites)
	RhoVector=array(0,c(2,MaxIts))
	
#the following includes member names in the code, however r cannot pass vectors of strings to fortran... need a fix...	
#ResultsRaw=.Fortran("FORMAT_PLINK_DATA",DataArray=as.integer(DataArray),SampleCoordinates=as.double(SampleCoordinates),MembersList=as.character(MembersList),PlinkFileName=as.character(paste(PlinkFileName,".ped",sep="")),LocationFileName=as.character(LocationFileName),SampleSites=as.integer(SampleSites),NumberSNPs=as.integer(NumberSNPs))

#ResultsRaw=.Fortran("FORMAT_PLINK_DATA",DataArray=as.integer(DataArray),SampleCoordinates=as.double(SampleCoordinates),PlinkFileName=as.character(paste(PlinkFileName,".ped",sep="")),LocationFileName=as.character(LocationFileName),SampleSites=as.integer(SampleSites),NumberSNPs=as.integer(NumberSNPs),PACKAGE="OriGen")
	#ResultsRaw=.Fortran("RTRIALS",DataArray=as.integer(DataArray),DataLength=as.integer(c(2,SampleSites,NumberSNPs)))
	
ResultsRaw=.Fortran("LEAVE_ONE_POP_OUT_CROSSVAL_SQUARE",PlinkFileName=as.character(PedFileName),LocationFileName=as.character(LocationFileName),NumberSNPs=as.integer(NumberSNPs),MaxIts=as.integer(MaxIts),MaxGridLength=as.integer(MaxGridLength),RhoVector=as.double(RhoVector),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),PACKAGE="OriGen")
	
#ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
#ResultsRaw$SampleCoordinates=aperm(array(ResultsRaw$SampleCoordinates,c(2,SampleSites)),c(2,1))
ResultsRaw$RhoVector=array(RhoVector,(c(2,MaxIts)))
	
	return(ResultsRaw)
}



FitAdmixedModelFindUnknowns<-function(DataArray,SampleCoordinates,UnknownData,MaxGridLength=20,RhoParameter=10,LambdaParameter=100.,MaskWater=TRUE){
#DataArray[Alleles,SampleSites,NumberSNPs] Gives the grouped data
#SampleCoordinates[SampleSites,2] gives the locations of the grouped data
#UnknownData[NumberUnknowns,NumberSNPs] gives the number of major alleles for the current unknown individual
#This function takes in the data, fits the model, and returns the allele frequency surfaces
	if(!.is.wholenumber(MaxGridLength)){
		stop("MaxGridLength must be an integer")
	}
	if(MaxGridLength<=1){
		stop("MaxGridLength must be greater than 1")
	}
	if(RhoParameter<=0){
		stop("RhoParameter must be greater than 0")
	}
	if(length(SampleCoordinates[1,])!=2){
		stop("SampleCoordinates should give the Long/Lat coordinates of the grouped data so it should only contain 2 columns")
	}
	NumberSNPs=length(DataArray[1,1,])
	SampleSites=length(DataArray[1,,1])
	GridLength=array(0,2)
	GridCoordinates=array(0.,dim=c(2,MaxGridLength))
	
	GridAndCoordResults=.Fortran("UPDATE_GRID_COORD_SQUARE2",GridCoordinates=as.double(GridCoordinates),SampleCoordinates=as.double(SampleCoordinates),GridLength=as.integer(GridLength),MaxGridLength=as.integer(MaxGridLength),SampleSites=as.integer(SampleSites),LambdaParameter=as.double(LambdaParameter),PACKAGE="OriGen")
	
	GridLength=GridAndCoordResults$GridLength
	GridCoordinates=array(GridAndCoordResults$GridCoordinates,c(2,MaxGridLength))
	print(GridLength)
	
	IsLand=array(TRUE,dim=c(GridLength[1],GridLength[2]))
	if(MaskWater){
		#change points on water to false here...
		IsLand=.LandArray(GridCoordinates,GridLength)
	}
	NumberUnknowns=length(UnknownData[,1])
	UnknownGrids=array(0,dim=c(GridLength[1],GridLength[2],NumberUnknowns))
	ResultsRaw=.Fortran("FITADMIXEDMODELFINDUNKNOWNS",AdmixtureFractions=as.double(UnknownGrids),DataArray=as.integer(DataArray),NumberSNPs=as.integer(NumberSNPs),GridLength=as.integer(GridLength),RhoParameter=as.double(RhoParameter),SampleSites=as.integer(SampleSites),MaxGridLength=as.integer(MaxGridLength),SampleCoordinates=as.double(SampleCoordinates),GridCoordinates=as.double(GridCoordinates),NumberUnknowns=as.integer(NumberUnknowns),UnknownData=as.integer(UnknownData),IsLand=as.logical(IsLand),PACKAGE="OriGen")
	
	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
	ResultsRaw$AdmixtureFractions=array(ResultsRaw$AdmixtureFractions,c(GridLength[1],GridLength[2],NumberUnknowns))
	ResultsRaw$SampleCoordinates=array(ResultsRaw$SampleCoordinates,c(SampleSites,2))
	ResultsRaw$UnknownData=array(ResultsRaw$UnknownData,c(NumberUnknowns,NumberSNPs))
	ResultsRaw$GridCoordinates=array(ResultsRaw$GridCoordinates,c(2,MaxGridLength))
	ResultsRaw$IsLand=array(ResultsRaw$IsLand,c(GridLength[1],GridLength[2]))
	
	return(ResultsRaw)
}




RankSNPsLRT<-function(DataArray){
#This function takes in the PED file along with a location file and outputs the Likelihood Ratio ranking
#of each SNP followed by the Likelihood Ratio statistic and the Informativeness for assignment by rosenberg et al..  Note that the statistic is compares the assumption
#that there is just a single global population vs several different sites.
	
	SampleSites=length(DataArray[1,,1])
	NumberSNPs=length(DataArray[1,1,])
	
	Rankings=1:NumberSNPs
	LRT=array(0,c(2,NumberSNPs))
	ResultsRaw=.Fortran("CALC_ALL_RANKINGS",DataArray=as.integer(DataArray),SampleSites=as.integer(SampleSites),NumberSNPs=as.integer(NumberSNPs),Rankings=as.integer(Rankings),LRT=as.double(LRT),PACKAGE="OriGen")
	
	ResultsRaw$DataArray=array(ResultsRaw$DataArray,c(2,SampleSites,NumberSNPs))
	ResultsRaw$LRT=array(ResultsRaw$LRT,c(2,NumberSNPs))
	return(ResultsRaw)
}




#-----------------------------------------------------------------------------------------------------

#The below functions will not be included in V1 of the R package....

#-----------------------------------------------------------------------------------------------------


#This function requires the maps package to work
.IsLand<-function(x.vec,y.vec){
	#require("maps")
	temp.vec=map.where(database="world",x.vec,y.vec)
	result.vec=x.vec*0+1
	for(i in 1:length(temp.vec)){
		if(is.na(temp.vec[i])){
			result.vec[i]=0
		}
	}
	return(result.vec)
}

#this function requires the maps package to work
.MaskWater<-function(GridCoordinates){
	#this short code checks whether the given coordinates are in water and outputs a matrix with 1 meaning land
	#and 0 meaning water...
	#GridCoordinates should be a matrix[x,2] where x is the number of grid points and the first 2 is Long,Lat
	ndiv=length(GridCoordinates[,1])
	latcount=0
	longcount=0
	for(i in 1:ndiv){
		if(GridCoordinates[i,2]>0.001){
			latcount=latcount+1
		}else if(GridCoordinates[i,2]< -0.001){
			latcount=latcount+1
		}
		if(GridCoordinates[i,1]>0.001){ 
			longcount=longcount+1
		}else if(GridCoordinates[i,1]< -0.001){
			longcount=longcount+1
		}
	}
	temp.mat=mat.or.vec(nc=longcount,nr=latcount)
	temp.mat[,]=1
	for(i in 1:longcount){
		temp.mat[,i]=.IsLand(rep(GridCoordinates[i,1],each=latcount),GridCoordinates[,2])
		}
	#write.table(temp.mat[latcount:1,],file="GridCoordSquare40Water.txt",append=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
	return(temp.mat)
}


#This function requires the maps package to work
.IsLandBool<-function(x.vec,y.vec){
	#require("maps")
	temp.vec=map.where(database="world",x.vec,y.vec)
	result.vec=x.vec
	result.vec[]=TRUE
	for(i in 1:length(temp.vec)){
		if(is.na(temp.vec[i])){
			result.vec[i]=FALSE
		}
	}
	return(result.vec)
}

.LandArray<-function(GridCoordinates,GridLength){
	#this short code checks whether the given coordinates are in water and outputs a matrix with 1 meaning land
	#and 0 meaning water...
	#GridCoordinates should be a matrix[2,x] where x is the number of grid points and the first 2 is Long,Lat
	#ndiv=length(GridCoordinates[1,])

	temp.mat=mat.or.vec(nr=GridLength[1],nc=GridLength[2])
	temp.mat[,]=TRUE
	for(i in 1:GridLength[1]){
		temp.mat[i,]=.IsLandBool(rep(GridCoordinates[1,i],each=GridLength[2]),GridCoordinates[2,1:GridLength[2]])
	}
	return(temp.mat)
}




#this function requires packages ggplot2 and maps to work.  Note that the vectors on the maps package is outdated particularly in europe
#An updated map can be downloaded from http://www.naturalearthdata.com/downloads/50m-cultural-vectors/


PlotAlleleFrequencySurface<-function(AlleleSurfaceOutput,SNPNumber=1,MaskWater=TRUE){
#GridCoordinates(2,MaxGridLength)
print("Note that the maps package used for vectors here is outdated, this is particularly true in Europe.") 
#require("maps")
#require("ggplot2")

TempHM=AlleleSurfaceOutput$AlleleFrequencySurfaces[SNPNumber,,]
for(i in 1:AlleleSurfaceOutput$GridLength[1]){
	TempHM[i,]=AlleleSurfaceOutput$GridCoordinates[1,i]
}
TempOb<-data.frame(Frequency=as.vector(AlleleSurfaceOutput$AlleleFrequencySurfaces[SNPNumber,,]),Long=as.vector(TempHM))
for(i in 1:AlleleSurfaceOutput$GridLength[2]){
	TempHM[,i]=AlleleSurfaceOutput$GridCoordinates[2,i]
}
TempOb$Lat=as.vector(TempHM)
TempOb$Land=.IsLand(TempOb$Long,TempOb$Lat)
subdata=subset(TempOb,Land==1)
#minp=min(subdata$Frequency)
minp=0
#maxp=max(subdata$Frequency)
maxp=1
if(MaskWater){
	subdata=subset(TempOb,Land==1)
	#minp=min(subdata$Frequency)
	#minp=0
	#maxp=max(subdata$Frequency)
	p<-ggplot(subset(TempOb,Land==1),aes(Long,Lat))
	} else {
	#minp=min(TempOb$Frequency)
	#minp=0
	#maxp=max(TempOb$Frequency)
	p<-ggplot(TempOb,aes(Long,Lat))
	}
p+	annotation_map(map_data("world"), fill=NA, colour = "white",asp=TRUE)+
	geom_tile(aes(fill=Frequency),colour=NA,alpha=1) +
	scale_fill_gradient(high = "#CFE8ED",low = "#0F4657",limits=c(minp,maxp)) +
	annotation_map(map_data("world",boundary=TRUE), fill=NA, colour = "black", bg=par(bg=NA)) + 
	ylab("Latitude") + ggtitle(paste0("Allele Frequency Surface SNP:",SNPNumber)) +
	xlab("Longitude")
}



PlotUnknownHeatMap<-function(HeatMapOutput,UnknownNumber=1,MaskWater=TRUE){
#GridCoordinates(2,MaxGridLength)
print("Note that the maps package used for vectors here is outdated, this is particularly true in Europe.") 
#require("maps")
#require("ggplot2")

TempHM=HeatMapOutput$UnknownGrids[,,UnknownNumber]
for(i in 1:HeatMapOutput$GridLength[1]){
	TempHM[i,]=HeatMapOutput$GridCoordinates[1,i]
}
TempOb<-data.frame(Probability=as.vector(HeatMapOutput$UnknownGrids[,,UnknownNumber]),Long=as.vector(TempHM))
for(i in 1:HeatMapOutput$GridLength[2]){
	TempHM[,i]=HeatMapOutput$GridCoordinates[2,i]
}
TempOb$Lat=as.vector(TempHM)
TempOb$Land=.IsLand(TempOb$Long,TempOb$Lat)
subdata=subset(TempOb,Land==1)
#minp=min(subdata$Probability)
minp=0
maxp=max(subdata$Probability)
if(MaskWater){
	subdata=subset(TempOb,Land==1)
	#minp=min(subdata$Probability)
	minp=0
	maxp=max(subdata$Probability)
	p<-ggplot(subset(TempOb,Land==1),aes(Long,Lat))
	} else {
	#minp=min(TempOb$Probability)
	minp=0
	maxp=max(TempOb$Probability)
	p<-ggplot(TempOb,aes(Long,Lat))
	}
p+	annotation_map(map_data("world"), fill=NA, colour = "white",asp=TRUE)+
	geom_tile(aes(fill=Probability),colour=NA,alpha=1) +
	scale_fill_gradient(high = "#CFE8ED",low = "#0F4657",limits=c(minp,maxp)) +
	annotation_map(map_data("world",boundary=TRUE), fill=NA, colour = "black", bg=par(bg=NA)) + 
	ylab("Latitude") + ggtitle(paste0("Heat Map Surface Individual:",UnknownNumber)) +
	xlab("Longitude")
}




PlotAdmixedSurface<-function(AdmixedOutput,UnknownNumber=1,MaskWater=TRUE){
#GridCoordinates(2,MaxGridLength)
print("Note that the maps package used for vectors here is outdated, this is particularly true in Europe.") 
#require("maps")
#require("ggplot2")

TempHM=AdmixedOutput$AdmixtureFractions[,,UnknownNumber]
#TempHM=AdmixedOutput$UnknownGrids[,,UnknownNumber]
for(i in 1:AdmixedOutput$GridLength[1]){
	TempHM[i,]=AdmixedOutput$GridCoordinates[1,i]
	}
TempOb<-data.frame(Fractions=as.vector(AdmixedOutput$AdmixtureFractions[,,UnknownNumber]),Long=as.vector(TempHM))
for(i in 1:AdmixedOutput$GridLength[2]){
	TempHM[,i]=AdmixedOutput$GridCoordinates[2,i]
	}
TempOb$Lat=as.vector(TempHM)
TempOb$Land=.IsLand(TempOb$Long,TempOb$Lat)
TempOb$Rounded<-round(TempOb$Fractions, digits=2)
subdata=subset(TempOb,Fractions>=0.01)

p<-ggplot(TempOb,aes(Long,Lat))
p+theme(panel.background = element_rect(fill = "lightskyblue1")) +
		annotation_map(map_data("world"), fill="darkolivegreen3", colour = "white",asp=TRUE)+
		annotation_map(map_data("world"),fill="NA",col="grey10") + 
		theme(legend.position="none") +
		geom_text(aes(label=Rounded),alpha=0,size=4) +
		geom_text(data=subdata,aes(Long,Lat,label=Rounded),alpha=1,size=4) +
		theme(legend.position="none") +
		ylab("Latitude")+xlab("Longitude")+ggtitle("Admixed Fractions")
}




