\name{missing.data}
\alias{missing.data}
\title{
Randomly input missing data points
}
\description{
Randomly replaces a set percentage of data points in a matrix of morphometric measurements with NA to simulate missing data. This is function RMD from Brown et al (2012). The amount of missing data can be chosen as an overall percentage of data (simple morphometric data) or specimens and can be constrained to a set of landmarks (for landmarks).
}
\usage{
missing.data(x, remperc, remsp = NULL, land.vec = NULL, land.identity = NULL)
}
\arguments{
  \item{x}{
A n X m matrix of morphometric data with n specimens and m variables. Or an array of geometric morphometrics landmarks (l X m X n)
}
  \item{remperc}{
The percentage of data to be removed from the matrix or array, expressed as a decimal (ex: 30 percent would be entered as 0.3)
}
  \item{remsp}{
The percentage of specimens to be removed from the array, expressed as a decimal (ex: 30 percent would be entered as 0.3)
}
  \item{land.vec}{
The number of landmarks to remove per specimen in an array. This can be a single value or vector with unique or repeating values.
}
  \item{land.identity}{
A vector to constrain the landmarks to chose from when assigning missing data. The values correspond to row numbers in an array.
}
}
\value{
Returns a n X m matrix or l X m X n array of morphometric data with missing variables input as NA
}
\references{
Brown, C., Arbour, J. and Jackson, D. 2012. Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology} 61(6):941-954.
}
\author{
J. Arbour and C. Brown
}
\seealso{
\code{\link{byclade}},\code{\link{obliterator}}
}
\examples{
data(dacrya)

#### remove 1 to 6 landmarks from 20% of specimens
dac.miss<-missing.data(dacrya,remsp=0.2,land.vec=c(1,2,3,4,5,6))
dac.miss
}

