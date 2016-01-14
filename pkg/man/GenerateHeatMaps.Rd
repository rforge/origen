\name{GenerateHeatMaps}
\alias{GenerateHeatMaps}
\title{
Fit OriGen microsatellite allele frequency surfaces
}
\description{
This function generates heat maps from OriGen microsatellite data output and then finds locations for unknown individuals..
}
\usage{
GenerateHeatMaps(FitModelOutput,UnknownDataArray,NumberLoci,MaskWater=TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
	
	\item{FitModelOutput}{This is the output from \code{\link{FitMultinomialModel}}.}
	
%	\item{PlinkFileName}{Base name of Plink PED file (i.e. without ".ped" or ".map") }

%  \item{LocationFileName}{Space or tab delimited text file with Longitude and Latitude coordinates for each individual listed in the 4th and 5th columns respectively.  Note that rows should correspond to the individuals in the Plink File.  Also, this file should have a header row.}
	
%	\item{DataArray}{An array giving the number of alleles grouped by sample sites for each SNP.  The dimension of this array is [MaxAlleles,SampleSites,NumberLoci].}
	
%	\item{SampleCoordinates}{This is an array which gives the longitude and latitude of each of the found sample sites.  The dimension of this array is [SampleSites,2], where the second dimension represents longitude and latitude respectively.}
	
	\item{UnknownDataArray}{This is an array which gives the alleles for the individuals of unknown origin.  The dimension of this array is [NumberUnknowns,2,NumberLoci], where 2 represents to 2 alleles each individual has at each locus.  Note that these should not be allele lengths but rather the allele number matching the dimension in DataArray.  Note that 0 or negative values here indicate unknown alleles and it is assumed that both are either known or unknown.}
	
%  \item{MaxGridLength}{An integer giving the maximum number of boxes to fill the longer side of the region.  Note that computation time increases quadratically as this number increases, but this number also should be high enough to separate different sample sites otherwise they will be binned together as a single site.}

	\item{NumberLoci}{This integer value gives the number of loci to include when generating the heat maps.  This is useful when generating heatmaps with multiple numbers of loci.}

%  \item{RhoParameter}{This is a real precision parameter weighting the amount of smoothing.  A higher value flattens out the surface while a lower value allows for more fluctuations.  The default value of 10 was used in our analysis and should prove a good starting point.  To choose a value by crossvalidation please see \code{\link{FindRhoParameterCrossValidation}}}
  
  \item{MaskWater}{If TRUE, this logical parameter restricts the heat maps to land areas only.}
 
}

\value{
List with the following components:
%\item{AlleleFrequencySurfaces}{An array giving the allele frequency for each allele, each coordinate, and each SNP.  The dimension of this array is [MaxAlleles,NumberLoci,NumberLongitudeDivisions,NumberLatitudeDivisions], where either NumberLongitudeDivisions or NumberLatitudeDivisions is equal to MaxGridLength(an input to this function) and the other is scaled so that the geodesic distance between points horizontally and vertically is equal.}

\item{AlleleFrequencySurfaces}{An array giving the allele frequency for each allele, each coordinate, and each SNP.  The dimension of this array is [MaxAlleles, NumberLoci, NumberLongitudeDivisions, NumberLatitudeDivisions], where either NumberLongitudeDivisions or NumberLatitudeDivisions is equal to MaxGridLength(an input to this function) and the other is scaled so that the geodesic distance between points horizontally and vertically is equal.}

\item{UnknownGrids}{An array giving the probability that an unknown individual comes from the given location.  The dimension of this array is [NumberLongitudeDivisions, NumberLatitudeDivisions, NumberUnknowns], where either NumberLongitudeDivisions or NumberLatitudeDivisions is equal to MaxGridLength(an input to this function) and the other is scaled so that the geodesic distance between points horizontally and vertically is equal.}

\item{DataArray}{An array giving the number alleles grouped by sample sites for each locus.  The dimension of this array is [MaxAlleles,SampleSites,NumberSNPs].}

\item{RhoParameter}{A real value showing the inputted RhoParameter value.}

\item{SampleSites}{This shows the integer number of sample sites found.}

\item{GridLength}{An array giving the number of longitudinal and latitudinal divisions.  The dimension of this array is [2], where the first number is longitude and the second is latitude.}

\item{MaxGridLength}{An integer giving the maximum number of boxes to fill the longer side of the region.  Note that computation time increases quadratically as this number increases, but this number also should be high enough to separate different sample sites otherwise they will be binned together as a single site.  This number was part of the inputs.}

\item{MaxAlleles}{This shows the maximum of AllelesAtLocus.  The maximum number of alleles at all loci.}

\item{NumberLoci}{This shows the integer number of loci found.}

\item{SampleCoordinates}{This is an array which gives the longitude and latitude of each of the found sample sites.  The dimension of this array is [SampleSites,2], where the second dimension represents longitude and latitude respectively.}

\item{GridCoordinates}{An array showing the corresponding coordinates for each longitude and latitude division.  The dimension of this array is [2,MaxGridLength], with longitude coordinates coming first and latitude second.  Note that one of these rows may not be filled entirely.  The associated output GridLength should be used to find the lengths of the two rows.  Rows not filled in entirely will contain zeroes at the end.}

\item{AllelesAtLocus}{This shows the integer vector of alleles found at each locus.}

\item{NumberUnknowns}{Integer number of unknown individuals found.}

\item{UnknownDataArray}{This is an array which gives the alleles for the individuals of unknown origin.  The dimension of this array is [NumberUnknowns,2,NumberLoci], where 2 represents to 2 alleles each individual has at each locus.  Note that these should not be allele lengths but rather the allele number matching the dimension in DataArray.}

%\item{UnknownPEDFile}{This shows the inputted PED file for the unknown individuals.}

%\item{NumberUnknowns}{This is an integer value showing the number of unknowns found in the UnknownPEDFile.}

%\item{UnknownData}{An array showing the unknown individuals genetic data.  The dimension of this array is [NumberUnknowns,NumberLoci].}

%\item{Membership}{This is an integer valued vector showing the group number of each member of the inputted known group.  The dimension of this array is [NumberKnown].}

%\item{NumberKnown}{This is an integer value showing the number of known found in the PlinkFileName.}

%\item{Rankings}{An integer valued vector giving the LRT based ranking of each SNP.  This can be used to reduce the number of SNPs to use for assignment if analysis takes too long.}

%\item{LRT}{This is a real valued array giving the Likelihood Ratio test statistic and the informativeness for assignment(Rosenberg) for each SNP.  The dimension of this array is [2,NumberLoci].}

}

\references{
Ranola J, Novembre J, Lange K (2014) Fast Spatial Ancestry via Flexible Allele Frequency Surfaces. Bioinformatics, in press.

}

\author{
John Michael Ranola, John Novembre, and Kenneth Lange
}

\seealso{
\code{\link{ConvertMicrosatData}} for converting Microsatellite data files into a format appropriate for analysis,	
	
\code{\link{ConvertPEDData}} for converting Plink PED files into a format appropriate for analysis,

\code{\link{FitMultinomialModel}} for fitting allele surfaces to the converted microsatellite data,

\code{\link{FitMultinomialModelFindUnknowns}} for fitting allele surfaces to the converted Microsatellite data,

\code{\link{PlotAlleleFrequencySurface}} for a quick way to plot the resulting allele frequency surfaces from \code{FitOriGenModel} or \code{GenerateHeatMaps},;

%\code{\link{ConvertUnknownPEDData}} for converting two Plink PED files (known and unknown)into a format appropriate for analysis,

%\code{\link{FitOriGenModelFindUnknowns}} for fitting allele surfaces to the converted data and finding the locations of the given unknown individuals,

%\code{\link{PlotUnknownHeatMap}} for a quick way to plot the resulting unknown heat map surfaces from \code{FitOriGenModelFindUnknowns},

%\code{\link{FitAdmixedFindUnknowns}} for fitting allele surfaces to the converted data and finding the locations of the given unknown individuals who may be admixed,

%\code{\link{PlotAdmixedSurface}} for a quick way to plot the resulting admixture surfaces from \code{FitAdmixedFindUnknowns},

%\code{\link{RankSNPsLRT}} for reducing the number of SNPs using a likelihood ratio test criteria or informativeness for assignment,

%\code{\link{FindRhoParamterCrossValidation}} for choosing an appropriate Rho parameter by way of crossvalidation,

}

\examples{

#this example not run because it takes longer than 5 secs
#note - type example(FunctionName, run.dontrun=TRUE) to run the example where FunctionName is
#the name of the function

\dontrun{

##Data generation
SampleSites=10
NumberLoci=4
MaxAlleles=4
NumberAllelesAtEachLocus=sample(2:MaxAlleles,NumberLoci,replace=TRUE)
TestData=array(0,dim=c(MaxAlleles,SampleSites,NumberLoci))
for(i in 1:NumberLoci){
	for(j in 1:NumberAllelesAtEachLocus[i]){
		TestData[j,,i]=sample(1:10,SampleSites,replace=TRUE)
	}
}
##This data is simulated in Europe which is around Longitude -9 to 38 and Latitude 34 to 60
TestCoordinates=array(0,dim=c(SampleSites,2))
TestCoordinates[,1]=runif(SampleSites,-9,38)
TestCoordinates[,2]=runif(SampleSites,34,60)

##This simulates the unknown data
NumberUnknowns=2
UnknownData=array(0,dim=c(NumberUnknowns,2,NumberLoci))
for(i in 1:NumberUnknowns){
	for(j in 1:NumberLoci){
		UnknownData[i,,j]=sample(1:NumberAllelesAtEachLocus[j],2)
	}
}

##MaxGridLength is the maximum number of boxes allowed 
##to span the region in either direction
##Note that this number was reduced to allow the example to run in less than 5 secs
##RhoParameter is a tuning constant
print("MaxGridLength is intentionally set really low for fast examples. 
	Meaningful results will most likely require a higher value.")

##Fits the allele frequency surfaces only
SurfaceTrials=FitMultinomialModel(TestData,TestCoordinates,
	MaxGridLength=20,RhoParameter=10)
str(SurfaceTrials)
##Plotting the model
PlotAlleleFrequencySurface(SurfaceTrials,LocusNumber=1,AlleleNumber=1,
	MaskWater=TRUE,Scale=FALSE)

##You can generate heatmaps of unknown individual's placements from with the allele
##surfaces using GenerateHeatMaps or use FitMultinomialModelFindUnknowns
HeatMapTrials=GenerateHeatMaps(SurfaceTrials,UnknownData,NumberLoci=NumberLoci)
##Plotting the unknown heat map
PlotUnknownHeatMap(HeatMapTrials,UnknownNumber=1,MaskWater=TRUE)
	
##Fitting the model and finding the unknown locations
#UnknownTrials=FitMultinomialModelFindUnknowns(TestData,TestCoordinates,
#	UnknownData,MaxGridLength=20,RhoParameter=10)
#str(UnknownTrials)
##Plotting the unknown heat map
#PlotUnknownHeatMap(UnknownTrials,UnknownNumber=1,MaskWater=TRUE)

##Fitting the admixed model
#AdmixedTrials=FitMultinomialAdmixedModelFindUnknowns(TestData,TestCoordinates,
#	UnknownData,MaxGridLength=10,RhoParameter=10)
##Plots the admixed surface disregarding fractions less than 0.01
#PlotAdmixedSurface(AdmixedTrials,UnknownNumber=1)

}

}


\keyword{localization}% __ONLY ONE__ keyword per line
\keyword{Ancestry}% __ONLY ONE__ keyword per line
