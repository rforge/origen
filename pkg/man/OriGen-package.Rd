\name{OriGen-package}
\alias{OriGen-package}
\alias{OriGen}
\docType{package}
\title{
Fast Spatial Ancestry via Flexible Allele Frequency Surfaces
}
\description{
This package primarily estimates allele frequency surfaces from point estimates.  It can also place individuals of unknown origin back onto the map with great accuracy.  Additionally, it can place admixed individuals by estimating contributing fractions at each location on a map. Lastly, it can rank SNPs by their ability to differentiate populations.
}
\details{
\tabular{ll}{
Package: \tab OriGen\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2013-10-13\cr
License: \tab GPL2 \cr
}
Index:
\itemize{
\item \code{\link{ConvertPEDData}}  This function converts Plink PED format files (PED/MAP) along with location files to the input required for OriGen.

\item \code{\link{ConvertUnknownPEDData}}    This function converts Plink PED format files (PED/MAP) along with location files to the input required for OriGen.  This differs from ConvertPEDData by its additional PED formatted input which contains the genotype information for unknown individuals.

\item \code{\link{FitOriGenModel}}    Fits the OriGen model and returns the allele frequency surfaces.  These surfaces can be plotted with the function PlotAlleleFrequencySurface.

\item \code{\link{PlotAlleleFrequencySurface}}    Plots a specified allele frequency surface from the output of FitOriGenModel.

\item \code{\link{FitOriGenModelFindUnknowns}}    Fits the OriGen model and places individuals of unknown origin onto the map.  This returns probability heat maps for each unknown individual.  These heat maps can be plotted with PlotUnknownHeatMap.

\item \code{\link{PlotUnknownHeatMap}}    Plots a specified unknown individuals heat map from the output of FitOriGenModelFindUnknowns.

\item \code{\link{FitAdmixedModelFindUnknowns}}    Fits the OriGen model and places unknown individuals who may be admixed onto the map.  Instead of returning a probability heat map for each individual, this returns admixture fractions at each location.  Note that many locations are 0.  This can be plotted with the function PlotAdmixedSurface.

\item \code{\link{PlotAdmixedSurface}}    Plots the admixture fractions of a specified individual from the output of FitAdmixedModelFindUnknowns.

\item \code{\link{RankSNPsLRT}}    This function takes a PED file along with a location file and outputs the likelihood ratio ranking of each SNP along with the LRT statistic and Rosenberg's informativeness for assignment.
}

}

\author{
John Michael Ranola, John Novembre, and Kenneth Lange

Maintainer: John Michael Ranola <ranolaj@uw.edu> 
}
\references{
Ranola J, Novembre J, and Lange K (2014) Fast Spatial Ancestry via Flexible Allele Frequency Surfaces.    Bioinformatics, in press.
}

\keyword{ Ancestry }
\keyword{ Admixture }
\keyword{ Localization }
