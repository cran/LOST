\name{est.reg}
\alias{est.reg}
\title{
A-priori size regression for missing data estimation
}
\description{
Estimates missing data using regression on a designated size variable. Any values of the size variable missing are estimated with the variable best correlated with size.
}
\usage{
est.reg(crocmiss, col_indep)
}
\arguments{
  \item{crocmiss}{
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
Brown, C., Arbour, J. and Jackson, D. (In Press) Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology}
}
\author{
J. Arbour and C. Brown
}
\seealso{
\code{\link{best.reg}}
}
