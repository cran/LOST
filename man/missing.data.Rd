\name{missing.data}
\alias{missing.data}
\title{
Randomly input missing data points
}
\description{
Randomly replaces a set percentage of data points in a matrix of morphometric measurements with NA to simulate missing data. This is function RMD from Brown et al. (In Press).
}
\usage{
missing.data(x, remperc)
}
\arguments{
  \item{x}{
A n X m matrix of morphometric data with n specimens and m variables
}
  \item{remperc}{
The percentage of data to be removed from the matrix, expressed as a decimal (ex: 30 percent would be entered as 0.3)
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
\code{\link{byclade}},\code{\link{obliterator}}
}
