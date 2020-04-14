\name{obliteratorGM}
\alias{obliteratorGM}
\title{
Simulate missing geometric morphometric landmarks with anatomical bias
}
\description{
This is the geometric morphometric implementation of the LOST function obliterator. This attempts to replicate specimens showing regional incompleteness. For each specimen, this function randomly selects one starting data point for removal. All subsequent data points have a probability of removal that is proportional to the inverse of the distance to all previously removed data points, based on the shape of that particular specimen (this differs from the linear morphometric implementation which requires a reference set of coordinates). For a complete mathematical description see Brown et al. (2012).
}
\usage{
obliteratorGM(x, remperc, expo=1)
}
\arguments{
  \item{x}{
A n X m matrix of morphometric data with n specimens and m variables. Or a l X 2 or 3 X n array of geometric morphometric coordinates, with l being the number of landmarks.
}
  \item{remperc}{
The percentage of data to be removed from the matrix, expressed as a decimal (ex: 30 percent would be entered as 0.3)
}
  \item{expo}{
An optional term for raising the denominator to an exponent, to increase or decrease the severity of the anatomical bias
}
}
\value{
Returns a n X m matrix of morphometric data with missing variables input as NA
}
\references{
Brown, C., Arbour, J. and Jackson, D. 2012. Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology} 61(6):941-954.
}
\author{
J. Arbour and C. Brown
}
\seealso{
\code{\link{missing.data}},\code{\link{byclade}},\code{\link{obliterator}}
}
