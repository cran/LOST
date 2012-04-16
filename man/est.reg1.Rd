\name{est.reg1}
\alias{est.reg1}
\title{
Internal function for best.reg
}
\description{
Internal function
}
\usage{
est.reg1(crocmiss, col_indep)
}
\arguments{
  \item{crocmiss}{
A n X m matrix of morphometric data with n specimens and m variables, containing some percentage of missing values input as NA
}
  \item{col_indep}{
The number of the column in which the independant size variable is stored. This column will be used to estimate missing values in the other columns.
}
}

\seealso{
\code{\link{best.reg}}
}
