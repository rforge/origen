\name{ConvertUnknownPEDData}
\alias{ConvertUnknownPEDData}
\title{
Plink PED file conversion for known and unknown data
}
\description{
This function converts two Plink PED/MAP files (one for the known samples and one with unknown locations) into the data format required for OriGen.
}
\usage{
ConvertUnknownPEDData(PlinkFileName,LocationFileName,PlinkUnknownFileName)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PlinkFileName}{Base name of Plink PED file (i.e. without ".ped" or ".map") containing the individuals with known locations.}

  \item{LocationFileName}{Space or tab delimited text file with Longitude and Latitude coordinates for each individual listed in the 4th and 5th columns respectively.  Note that rows should correspond to the individuals in the Plink File.  Also, this file should have a header row.}
  
  \item{PlinkUnknownFileName}{Base name of Plink PED file (i.e. without ".ped" or ".map") containing the individuals with unknown locations.}
  
 
}

\value{
List with the following components:

\item{DataArray}{An array giving the number of major/minor SNPs (defined as the most occuring in the dataset) grouped by sample sites for each SNP.  The dimension of this array is [2,SampleSites,NumberSNPs].}

\item{SampleCoordinates}{This is an array which gives the longitude and latitude of each of the found sample sites.  The dimension of this array is [SampleSites,2], where the second dimension represents longitude and latitude respectively.}

\item{PlinkFileName}{This shows the inputted PlinkFileName with ".ped" attached.}

\item{LocationFile}{This shows the inputted LocationFileName.}

\item{SampleSites}{This shows the integer number of sample sites found.}

\item{NumberSNPs}{This shows the integer number of SNPs found.}

\item{UnknownPEDFile}{This shows the inputted PED file for the unknown individuals.}

\item{NumberUnknowns}{This is an integer value showing the number of unknowns found in the UnknownPEDFile.}

\item{UnknownData}{An array showing the unknown individuals genetic data.  The dimension of this array is [NumberUnknowns,NumberSNPs].}

\item{Membership}{This is an integer valued vector showing the group number of each member of the inputted known group.  The dimension of this array is [NumberKnown].}

\item{NumberKnown}{This is an integer value showing the number of known found in the PlinkFileName.}
}

\references{
Ranola J, Novembre J, Lange K (2014) Fast Spatial Ancestry via Flexible Allele Frequency Surfaces. Bioinformatics, in press.

}

\author{
John Michael Ranola, John Novembre, and Kenneth Lange
}

\seealso{

%\code{\link{ConvertPEDData}} for converting Plink PED files into a format appropriate for analysis,

%\code{\link{FitOriGenModel}} for fitting allele surfaces to the converted data,

%\code{\link{PlotAlleleFrequencySurface}} for a quick way to plot the resulting allele frequency surfaces from \code{FitOriGenModel},

\code{\link{ConvertUnknownPEDData}} for converting two Plink PED files (known and unknown)into a format appropriate for analysis,

\code{\link{FitOriGenModelFindUnknowns}} for fitting allele surfaces to the converted data and finding the locations of the given unknown individuals,

\code{\link{PlotUnknownHeatMap}} for a quick way to plot the resulting unknown heat map surfaces from \code{FitOriGenModelFindUnknowns},;

%\code{\link{FitAdmixedFindUnknowns}} for fitting allele surfaces to the converted data and finding the locations of the given unknown individuals who may be admixed,

%\code{\link{PlotAdmixedSurface}} for a quick way to plot the resulting admixture surfaces from \code{FitAdmixedFindUnknowns},

%\code{\link{RankSNPsLRT}} for reducing the number of SNPs using a likelihood ratio test criteria or informativeness for assignment,

%\code{\link{FindRhoParamterCrossValidation}} for choosing an appropriate Rho parameter by way of crossvalidation,

}

\examples{

#Note that Plink files "10SNPs.ped", "10SNPs.map" and also "Locations.txt" 
#are included in the data folder of the OriGen package with ".txt" appended to the Plink files.
#Please remove ".txt" and navigate to the appropriate location   
#before testing the following commands.
#Note that this was done to allow inclusion of the test data in the package.

\dontrun{trials3=ConvertUnknownPEDData("10SNPs","Locations.txt",""10SNPs"")}
\dontrun{str(trials3)}
MaxGridLength=30
RhoParameter=10
\dontrun{trials4=FitOriGenModelFindUnknowns(trials3$DataArray,trials3$SampleCoordinates,
	trials3$UnknownData[1:2,],MaxGridLength,RhoParameter)}
\dontrun{PlotUnknownHeatMap(trials4,UnknownNumber=1,MaskWater=TRUE)}

}


\keyword{Conversion}% __ONLY ONE__ keyword per line
\keyword{Plink}% __ONLY ONE__ keyword per line
\keyword{Files}% __ONLY ONE__ keyword per line
\keyword{PED}% __ONLY ONE__ keyword per line
