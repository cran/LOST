\name{MissingGeoMorph}
\alias{MissingGeoMorph}
\title{
Estimate missing landmark data
}
\description{
This function provides several options for estimating landmark data (details of which can be found in the references below). The function first aligned the landmarks using Procrustes superimposition (\code{\link{align.missing}}).
}
\usage{
MissingGeoMorph(x, nlandmarks, method = "BPCA")
}
\arguments{
  \item{x}{
A n* l X 2 matrix of coordinate data, where n is the number of specimens and l is the number of landmarks. All landmarks from one specimen should be grouped together. Missing values should be given as NA
}
  \item{nlandmarks}{
The number of landmarks per specimen.
}
  \item{method}{
Four methods are provided for estimating missing landmark data: 1) "BPCA" - Bayesian principal component analysis, 2) "mean" - mean substitution, 3) "reg" - values are estimated based on the highest correlated variable available, and 4) "TPS" - thin plate spline interpolation. See Arbour and Brown (In Press) for a comparison of the performance of each of these methods.
}
}

\value{
Returns an n * l X 2 matrix of coordinate data, with missing values imputed. Landmarks have been aligned and are given in the original shape space.
}
\references{
Arbour, J. and Brown, C. \emph{In Press}. Incomplete specimens in Geometric Morphometric Analyses. \emph{Methods in Ecology and Evolution}

Brown, C., Arbour, J. and Jackson, D. 2012. Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology} 61(6):941-954.
}
\author{
J. Arbour
}

\seealso{
\code{\link{align.missing}}, \code{\link{missing.specimens}}
}
\examples{

data(dacrya)

#### make some specimens incomplete
dac.miss<-missing.specimens(dacrya,10,c(1,2,3,4,5,6),16)
dac.miss

####estimate missing landmarks using mean substituion
dac.new<-MissingGeoMorph(dac.miss,16,method="mean")
dac.new

}