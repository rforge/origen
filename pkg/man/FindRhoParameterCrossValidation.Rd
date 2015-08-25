\name{FindRhoParameterCrossValidation}
\alias{FindRhoParameterCrossValidation}
\title{
Finds the appropriate value of the Rho parameter via crossvalidation.
}
\description{
This function finds the appropriate value of the tuning constant, RhoParameter, via a leave one sample site out cross validation.
}
\usage{
FindRhoParameterCrossValidation(PlinkFileName,LocationFileName,MaxIts=6,MaxGridLength=20)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
	%\item{AlleleSurfaceOutput}{The output of \code{\link{FitOriGenModel}}}
	%\item{SNPNumber}{Integer indicating the SNP number to plot.}
	%\item{MaskWater}{Logical value that if true removes water from the plotted regions.}
	
	\item{PlinkFileName}{Base name of Plink PED file (i.e. without ".ped" or ".map") }

  \item{LocationFileName}{Space or tab delimited text file with Longitude and Latitude coordinates for each individual listed in the 4th and 5th columns respectively.  Note that rows should correspond to the individuals in the Plink File.  Also, this file should have a header row.}
  
  \item{MaxIts}{An integer giving the number of iterations before selecting the rho parameter. Note that this is a long process so it is best to start small.}
  	
	%\item{DataArray}{An array giving the number of major/minor SNPs (defined as the most occuring in the dataset) grouped by sample sites for each SNP.  The dimension of this array is [2,SampleSites,NumberSNPs].}
	
	%\item{SampleCoordinates}{This is an array which gives the longitude and latitude of each of the found sample sites.  The dimension of this array is [SampleSites,2], where the second dimension represents longitude and latitude respectively.}
	
	%\item{UnknownData}{An array showing the unknown individuals genetic data.  The dimension of this array is [NumberUnknowns,NumberSNPs].}
	
  \item{MaxGridLength}{An integer giving the maximum number of boxes to fill the longer side of the region.  Note that computation time increases quadratically as this number increases, but this number also should be high enough to separate different sample sites otherwise they will be binned together as a single site.}

  %\item{RhoParameter}{This is a real precision parameter weighting the amount of smoothing.  A higher value flattens out the surface while a lower value allows for more fluctuations.  The default value of 10 was used in our analysis and should prove a good starting point.  To choose a value by crossvalidation please see \code{\link{FindRhoParameterCrossValidation}}}
 
}

\value{
	%This outputs a plot of the allele frequency surface on a map.
List with the following components:
%\item{AlleleFrequencySurfaces}{An array giving the allele frequency for each coordinate and each SNP.  The dimension of this array is [NumberSNPs,NumberLongitudeDivisions,NumberLatitudeDivisions], where either NumberLongitudeDivisions or NumberLatitudeDivisions is equal to MaxGridLength(an input to this function) and the other is scaled so that the geodesic distance between points horizontally and vertically is equal.}

%\item{UnknownGrids}{An array giving the probability that an unknown individual comes from the given location.  The dimension of this array is [NumberLongitudeDivisions,NumberLatitudeDivisions,NumberUnknowns], where either NumberLongitudeDivisions or NumberLatitudeDivisions is equal to MaxGridLength(an input to this function) and the other is scaled so that the geodesic distance between points horizontally and vertically is equal.}

%\item{AdmixtureFractions}{An array giving the admixture fraction from the given location.  In other words this is the fractional contribution of the location to the unknown individuals genetic data.  The dimension of this array is [NumberLongitudeDivisions,NumberLatitudeDivisions,NumberUnknowns], where either NumberLongitudeDivisions or NumberLatitudeDivisions is equal to MaxGridLength(an input to this function) and the other is scaled so that the geodesic distance between points horizontally and vertically is equal.}

\item{PlinkFileName}{This shows the inputted PlinkFileName with ".ped" attached.}

\item{LocationFile}{This shows the inputted LocationFileName.}

\item{NumberSNPs}{This shows the integer number of SNPs found.}

\item{MaxIts}{An integer giving the number of iterations before selecting the rho parameter. Note that this is a long process so it is best to start small.  This number is inputted into the function.}

\item{MaxGridLength}{An integer giving the maximum number of boxes to fill the longer side of the region.  Note that computation time increases quadratically as this number increases, but this number also should be high enough to separate different sample sites otherwise they will be binned together as a single site.  This number was part of the inputs.}

%\item{DataArray}{An array giving the number of major/minor SNPs (defined as the most occuring in the dataset) grouped by sample sites for each SNP.  The dimension of this array is [2,SampleSites,NumberSNPs].}

\item{RhoVector}{An array giving the tested values of RhoParameter along with the resulting cross validation results where lower is better.}

\item{GridLength}{An array giving the number of longitudinal and latitudinal divisions.  The dimension of this array is [2], where the first number is longitude and the second is latitude.}

\item{RhoParameter}{A real value showing the best RhoParameter value found.}

\item{SampleSites}{This shows the integer number of sample sites found.}

%\item{MaxGridLength}{An integer giving the maximum number of boxes to fill the longer side of the region.  Note that computation time increases quadratically as this number increases, but this number also should be high enough to separate different sample sites otherwise they will be binned together as a single site.  This number was part of the inputs.}

%\item{SampleCoordinates}{This is an array which gives the longitude and latitude of each of the found sample sites.  The dimension of this array is [SampleSites,2], where the second dimension represents longitude and latitude respectively.}

%\item{GridCoordinates}{An array showing the corresponding coordinates for each longitude and latitude division.  The dimension of this array is [2,MaxGridLength], with longitude coordinates coming first and latitude second.  Note that one of these rows may not be filled entirely.  The associated output GridLength should be used to find the lengths of the two rows.  Rows not filled in entirely will contain zeroes at the end.}

%\item{UnknownPEDFile}{This shows the inputted PED file for the unknown individuals.}

%\item{NumberUnknowns}{This is an integer value showing the number of unknowns found in the UnknownPEDFile.}

%\item{UnknownData}{An array showing the unknown individuals genetic data.  The dimension of this array is [NumberUnknowns,NumberSNPs].}

%\item{Membership}{This is an integer valued vector showing the group number of each member of the inputted known group.  The dimension of this array is [NumberKnown].}

%\item{NumberKnown}{This is an integer value showing the number of known found in the PlinkFileName.}

%\item{Rankings}{An integer valued vector giving the LRT based ranking of each SNP.  This can be used to reduce the number of SNPs to use for assignment if analysis takes too long.}

%\item{LRT}{This is a real valued array giving the Likelihood Ratio test statistic and the informativeness for assignment(Rosenberg) for each SNP.  The dimension of this array is [2,NumberSNPs].}

%\item{IsLand}{This is a logical valued array that is TRUE when the given coordinates are over land and FALSE when over water.  The dimension of this array is [GridLength[1],GridLength[2]].}

}

\references{
Ranola J, Novembre J, Lange K (2014) Fast Spatial Ancestry via Flexible Allele Frequency Surfaces. Bioinformatics, in press.

}

\author{
John Michael Ranola, John Novembre, and Kenneth Lange
}

\seealso{
\code{\link{ConvertPEDData}} for converting Plink PED files into a format appropriate for analysis,

\code{\link{FitOriGenModel}} for fitting allele surfaces to the converted data,

\code{\link{PlotAlleleFrequencySurface}} for a quick way to plot the resulting allele frequency surfaces from \code{FitOriGenModel},

\code{\link{ConvertUnknownPEDData}} for converting two Plink PED files (known and unknown)into a format appropriate for analysis,

\code{\link{FitOriGenModelFindUnknowns}} for fitting allele surfaces to the converted data and finding the locations of the given unknown individuals,

\code{\link{PlotUnknownHeatMap}} for a quick way to plot the resulting unknown heat map surfaces from \code{FitOriGenModelFindUnknowns},;

\code{\link{FitAdmixedModelFindUnknowns}} for fitting allele surfaces to the converted data and finding the locations of the given unknown individuals who may be admixed,

\code{\link{PlotAdmixedSurface}} for a quick way to plot the resulting admixture surfaces from \code{FitAdmixedFindUnknowns},

\code{\link{RankSNPsLRT}} for reducing the number of SNPs using a likelihood ratio test criteria or informativeness for assignment,

%\code{\link{FindRhoParamterCrossValidation}} for choosing an appropriate Rho parameter by way of crossvalidation,

}

\examples{
#Note that Plink files "10SNPs.ped", "10SNPs.map" and also "Locations.txt" 
#are included in the data folder of the OriGen package.  
#Please navigate to the appropriate location before testing 
#the following commands.

\dontrun{trials5=FindRhoParameterCrossValidation("10SNPs","Locations.txt",
	MaxIts=4,MaxGridLength=20)}
\dontrun{trials5}

}

\keyword{Crossvalidation}% __ONLY ONE__ keyword per line
\keyword{Tuning}% __ONLY ONE__ keyword per line
