\name{est.reg}
\alias{est.reg}
\title{
A-priori size regression for missing data estimation
}
\description{
Estimates missing data using regression on a designated size variable. Any values of the size variable missing are estimated with the variable best correlated with size.
}
\usage{
est.reg(x, col_indep)
}
\arguments{
  \item{x}{
A n X m matrix of morphometric data with n specimens and m variables, containing some percentage of missing values input as NA
}
  \item{col_indep}{
The number of the column in which the independant size variable is stored. This column will be used to estimate missing values in the other columns.
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
\code{\link{best.reg}}
}
\examples{

data(crocs)

## remove 30% of data points
croc.miss<-missing.data(crocs,0.3)
croc.miss

## assume col 1 is the size variable
croc.new<-est.reg(croc.miss,1)
croc.new
}

