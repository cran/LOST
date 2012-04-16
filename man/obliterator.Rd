\name{obliterator}
\alias{obliterator}
\title{
Simulate missing morphometric data with anatomical bias
}
\description{
This function simulates the effect of proximity between measurements in morphometric data on the distribution of missing values. This attempts to replicate specimens showing regional deformation or incompleteness. From a morphometric dataset, this function selects a number of specimens to have data points removed from and a number of measurements to remove from each of these specimens based on the distribution of missing data produced by \code{\link{missing.data}}. For each specimen, this function randomly selects one starting data point for removal. All subsequent data points have a probability of removal that is proportional to the inverse of the distance to all previously removed data points, based on a reference set of landmarks (matrix 'distances'). For a complete mathematical description see Brown et al. (In Press).
}
\usage{
obliterator(x, remperc, distances, expo=1)
}
\arguments{
  \item{x}{
A n X m matrix of morphometric data with n specimens and m variables
}
  \item{remperc}{
The percentage of data to be removed from the matrix, expressed as a decimal (ex: 30 percent would be entered as 0.3)
}
  \item{distances}{
An m X 6 matrix include the start and end points (landmarks) for each morphometric measurement from a reference specimen (3D). The data in each row is ordered as x1,x2,y1,y2,z1,z2.
}
  \item{expo}{
An optional term for raising the denominator to an exponent, to increase or decrease the severity of the anatomical bias
}
}
\value{
Returns a n X m matrix of morphometric data with missing variables input as NA
}
\references{
Brown, C., Arbour, J. and Jackson, D. (In Press) Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology}
}
\author{
J. Arbour and C. Brown
}
\seealso{
\code{\link{missing.data}},\code{\link{byclade}}
}
