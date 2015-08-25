\name{ConvertPEDData}
\alias{ConvertPEDData}
\title{
Plink PED file conversion
}
\description{
This function converts a Plink PED/MAP file into the data format required for OriGen.
}
\usage{
ConvertPEDData(PlinkFileName,LocationFileName)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PlinkFileName}{Base name of Plink PED file (i.e. without ".ped" or ".map") 
}

  \item{LocationFileName}{Space or tab delimited text file with Longitude and Latitude coordinates for each individual listed in the 4th and 5th columns respectively.  Note that rows should correspond to the individuals in the Plink File.  Also, this file should have a header row.}
 
}

\value{
List with the following components:

\item{DataArray}{An array giving the number of major/minor SNPs (defined as the most occuring in the dataset) grouped by sample sites for each SNP.  The dimension of this array is [2,SampleSites,NumberSNPs].}

\item{SampleCoordinates}{This is an array which gives the longitude and latitude of each of the found sample sites.  The dimension of this array is [SampleSites,2], where the second dimension represents longitude and latitude respectively.}

\item{PlinkFileName}{This shows the inputted PlinkFileName with ".ped" attached.}

\item{LocationFile}{This shows the inputted LocationFileName.}

\item{SampleSites}{This shows the integer number of sample sites found.}

\item{NumberSNPs}{This shows the integer number of SNPs found.}


}

\references{
Ranola J, Novembre J, Lange K (2014) Fast Spatial Ancestry via Flexible Allele Frequency Surfaces. Bioinformatics, in press.

}

\author{
John Michael Ranola, John Novembre, and Kenneth Lange
}

\seealso{

\code{\link{FitOriGenModel}} for fitting allele surfaces to the converted data,

\code{\link{PlotAlleleFrequencySurface}} for a quick way to plot the resulting allele frequency surfaces from \code{FitOriGenModel}

\code{\link{ConvertUnknownPEDData}} for converting a known and unknown PED files (2 separate files) into the format required for OriGen (Note that this is what you want if you want to place unknown individuals back on the map);

}

\examples{

#Note that Plink files "10SNPs.ped", "10SNPs.map" and also "Locations.txt" 
#are included in the data folder of the OriGen package with ".txt" appended to the Plink files.
#Please remove ".txt" and navigate to the appropriate location   
#before testing the following commands.
#Note that this was done to allow inclusion of the test data in the package.

\dontrun{trials=ConvertPEDData("10SNPs","Locations.txt")}
\dontrun{str(trials)}
MaxGridLength=20
RhoParameter=10
\dontrun{trials2=FitOriGenModel(trials$DataArray,trials$SampleCoordinates,
	MaxGridLength,RhoParameter)}
\dontrun{PlotAlleleFrequencySurface(trials2)}

}


\keyword{Conversion}% __ONLY ONE__ keyword per line
\keyword{Plink}% __ONLY ONE__ keyword per line
\keyword{Files}% __ONLY ONE__ keyword per line
\keyword{PED}% __ONLY ONE__ keyword per line

