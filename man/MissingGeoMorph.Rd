\name{MissingGeoMorph}
\alias{MissingGeoMorph}
\title{
Estimate missing landmark data
}
\description{
This function provides several options for estimating landmark data (details of which can be found in the references below). The function first alignes the landmarks using Procrustes superimposition (\code{\link{align.missing}}). Both 2D and 3D coordinates can be accommodated.
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
Four methods are provided for estimating missing landmark data: 1) "BPCA" - Bayesian principal component analysis, 2) "mean" - mean substitution, 3) "reg" - values are estimated based on the most strongly correlated variable available, and 4) "TPS" - thin plate spline interpolation (only available for 2D). See Arbour and Brown (2014) for a comparison of the performance of each of these methods.
}
}

\value{
Returns an n * l X 2 (or 3) matrix of coordinate data, with missing values imputed. Landmarks have been aligned and are given in the original shape space.
}
\references{
Arbour, J. and Brown, C. 2014. Incomplete specimens in Geometric Morphometric Analyses. \emph{Methods in Ecology and Evolution} 5(1):16-26.

Brown, C., Arbour, J. and Jackson, D. 2012. Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology} 61(6):941-954.
}
\author{
J. Arbour
}

\seealso{
\code{\link{align.missing}}, \code{\link{missing.specimens}}
}
