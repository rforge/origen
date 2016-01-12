\name{CalcFractionsMultiLoglik}
\alias{CalcFractionsMultiLoglik}
\title{
Calculates the loglikelihood for placing a sample 100 percent back into its own sample site
}
\description{
This function takes the UnknownDataArray which contains allelelic information
	for individuals WITHIN a single sample site and calculates the resulting
	fraction loglikelihood for placing all individuals 100 percent back into their site
}
\usage{
CalcFractionsMultiLoglik(UnknownDataArray,LambdaParameter=100)

}
%- maybe also 'usage' for other objects documented here.
\arguments{

	\item{UnknownDataArray}{An array showing the unknown individuals genetic data.  It lists the two allele numbers of the unknown data.  The dimension of this array is [NumberUnknowns,2,NumberLoci].}

  
	\item{LambdaParameter}{This is a real precision parameter weighting the admixture fractions algorithm.  For the most part, this does not need to be changed as it seems to only affect the time to convergence.  Default is 100.}
  
 
}

\value{

An array giving the penalized loglikelihood resulting from placing each unknown individual 100 percent back into his own sample site.  The length of this array is [NumberUnknowns].



}

\references{
Ranola J, Novembre J, Lange K (2014) Fast Spatial Ancestry via Flexible Allele Frequency Surfaces. Bioinformatics, in press.

}

\author{
John Michael Ranola, John Novembre, and Kenneth Lange
}

\seealso{
\code{\link{FitMultinomialAdmixedModelFindUnknowns}} for getting loglikelihoods of unknown individuals placed into chosen regions.

%\code{\link{ConvertPEDData}} for converting Plink PED files into a format appropriate for analysis,

%\code{\link{FitOriGenModel}} for fitting allele surfaces to the converted data,

%\code{\link{PlotAlleleFrequencySurface}} for a quick way to plot the resulting allele frequency surfaces from \code{FitOriGenModel},

%\code{\link{ConvertUnknownPEDData}} for converting two Plink PED files (known and unknown)into a format appropriate for analysis,

%\code{\link{FitOriGenModelFindUnknowns}} for fitting allele surfaces to the converted data and finding the locations of the given unknown individuals,

%\code{\link{PlotUnknownHeatMap}} for a quick way to plot the resulting unknown heat map surfaces from \code{FitOriGenModelFindUnknowns},;

%\code{\link{CalcFractionsMultiLoglik}} for fitting allele surfaces to the converted data and finding the locations of the given unknown individuals who may be admixed,

%\code{\link{PlotAdmixedSurface}} for a quick way to plot the resulting admixture surfaces from \code{FitAdmixedFindUnknowns},

%\code{\link{RankSNPsLRT}} for reducing the number of SNPs using a likelihood ratio test criteria or informativeness for assignment,

%\code{\link{FindRhoParamterCrossValidation}} for choosing an appropriate Rho parameter by way of crossvalidation,

}

\examples{

#Data generation
NumberUnknowns = 50
NumberLoci = 10
TestUnknownDataArray=array(sample(1:5,2*NumberUnknowns*NumberLoci,replace=TRUE)
	,dim=c(NumberUnknowns,2,NumberLoci))

CalcFractionsMultiLoglik(TestUnknownDataArray)



}


\keyword{Localization}% __ONLY ONE__ keyword per line
\keyword{Ancestry}% __ONLY ONE__ keyword per line
\keyword{Admixture}% __ONLY ONE__ keyword per line
