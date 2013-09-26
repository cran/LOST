\name{best.reg}
\alias{best.reg}
\title{
Estimate missing morphometric data with a highly correlated variable
}
\description{
Estimates missing morphometric using regression on the most highly correlated morphological variable available 
}
\usage{
best.reg(x)
}
\arguments{
  \item{x}{
A n X m matrix of morphometric data with n specimens and m variables, containing some percentage of missing values input as NA
}
}
\value{
Returns a n X m matrix containing both the original morphometric values as well as estimates for all previously missing values.
}
\references{
Brown, C., Arbour, J. and Jackson, D. 2012. Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology} 61(6):941-954.
}
\author{
J. Arbour and C. Brown
}
\seealso{
\code{\link{est.reg}}
}
\examples{

data(crocs)

### remove 30% of data points
croc.miss<-missing.data(crocs,0.3)
croc.miss

### estimating missing points
croc.new<-best.reg(croc.miss)
croc.new
}