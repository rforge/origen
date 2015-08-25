\name{ConvertMicrosatData}
\alias{ConvertMicrosatData}
\title{
Microsatellite file conversion for known and unknown data
}
\description{
This function converts two Microsatellite data files (one for the genotypes and one for locations) into the data format required for OriGen.
}
\usage{
ConvertMicrosatData(DataFileName,LocationFileName)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataFileName}{Name of file containing the genotypes of the various locations.  The columns here would be LocationName, LocationNumber, Locus1, Locus2, etc.  Each individual would take up 2 rows (one for each allele) with the same LocationName and LocationNumber.  The value under Locus would be the length of the allele of that individual.  Note that unknown individuals should have location number "-1".}

  \item{LocationFileName}{Space or tab delimited text file with the location information for the individuals.  The columns are LocationName, LocationNumber, Latitude, and Longitude.  Note that the first two columns must be in the same order as the FileName.}
 
}

\value{
List with the following components:

\item{DataArray}{An array giving the number alleles grouped by sample sites for each locus.  The dimension of this array is [MaxAlleles,SampleSites,NumberSNPs].}

\item{SampleCoordinates}{This is an array which gives the longitude and latitude of each of the found sample sites.  The dimension of this array is [SampleSites,2], where the second dimension represents longitude and latitude respectively.}

\item{AllelesAtLocus}{This shows the integer vector of alleles found at each locus.}

\item{MaxAlleles}{This shows the maximum of AllelesAtLocus.  The maximum number of alleles at all loci.}

\item{SampleSites}{This shows the integer number of sample sites found.}

\item{NumberLoci}{This shows the integer number of loci found.}

\item{NumberUnknowns}{This is an integer value showing the number of unknowns found.}

\item{UnknownDataArray}{An array showing the unknown individuals genetic data.  The dimension of this array is [NumberUnknowns,2,NumberLoci].}

\item{LocationNames}{This is a list of all the LocationNames (The first column of the input files).}

\item{DataFileName}{This shows the inputted DataFileName.}

\item{LocationFileName}{This shows the inputted LocationFileName.}
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

\code{\link{ConvertMicrosatData}} for converting Microsatellite data files into a format appropriate for analysis,	
	
\code{\link{ConvertPEDData}} for converting Plink PED files into a format appropriate for analysis,

\code{\link{FitMultinomialModel}} for fitting allele surfaces to the converted Microsatellite data,

\code{\link{PlotAlleleFrequencySurface}} for a quick way to plot the resulting allele frequency surfaces from \code{FitOriGenModel} or \code{FitMultinomialModel},;

%\code{\link{FitAdmixedFindUnknowns}} for fitting allele surfaces to the converted data and finding the locations of the given unknown individuals who may be admixed,

%\code{\link{PlotAdmixedSurface}} for a quick way to plot the resulting admixture surfaces from \code{FitAdmixedFindUnknowns},

%\code{\link{RankSNPsLRT}} for reducing the number of SNPs using a likelihood ratio test criteria or informativeness for assignment,

%\code{\link{FindRhoParamterCrossValidation}} for choosing an appropriate Rho parameter by way of crossvalidation,

}

\examples{

#Note that sample files MicrosatTrialDataSmall.txt and 
#LocationTrialDataSmall.txt are included in data for formatting.
#Note that this was done to allow inclusion of the test data in the package.

\dontrun{MicrosatDataSmall=ConvertMicrosatData("MicrosatTrialDataSmall.txt",
		"LocationTrialDataSmall.txt")}
\dontrun{str(MicrosatDataSmall)}
\dontrun{MicrosatAnalysisSmall=FitMultinomialModel(MicrosatDataSmall$DataArray,
		MicrosatDataSmall$SampleCoordinates,MaxGridLength=20)}
\dontrun{str(MicrosatAnalysisSmall)}
\dontrun{PlotAlleleFrequencySurface(MicrosatAnalysisSmall)}


}


\keyword{Conversion}% __ONLY ONE__ keyword per line
\keyword{Plink}% __ONLY ONE__ keyword per line
\keyword{Files}% __ONLY ONE__ keyword per line
\keyword{PED}% __ONLY ONE__ keyword per line
