\name{how.many.missing}
\alias{how.many.missing}
\title{
Calculate the percentage of missing morphometric data
}
\description{
Calculates the percentage of morphometric data points that have been replaced with 'NA' by functions such as \code{\link{missing.data}}, \code{\link{byclade}} or \code{\link{obliterator}} from LOST. Used to verify the amount of missing data inputted into complete morphometric matrices.
}
\usage{
how.many.missing(x)
}
\arguments{
  \item{x}{
A n X m matrix of morphometric data with n specimens and m variables containing some percentage of missing data
}
}
\value{
Returns the percentage (as a decimal) of missing data points present in x
}
\references{
Brown, C., Arbour, J. and Jackson, D. (In Press) Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology}
}
\author{
J. Arbour and C. Brown
}
\seealso{
\code{\link{missing.data}}
}
