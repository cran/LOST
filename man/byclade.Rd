\name{byclade}
\alias{byclade}
\title{
Simulate missing morphometric data with taxonomic bias
}
\description{
This function simulates higher frequency of missing data points in groups that are less numerically well represented in the whole sample, relative to other group. These groups may represent taxa (as used in Brown et al., In Press), but may also represent any other group of interest (e.g. populations, trials, subsamples, etc.). From a morphometric dataset, this function selects a number of specimens to have data points removed from and a number of measurements to remove from each of these specimens based on the distribution of missing data produced by \code{\link{missing.data}}. A vector containing the number of measurements to remove from each specimen is produced and sorted into descending order. Specimens are then sampled without replacement with a probability relative to the sum of the entire sample sizes divided by the number of specimens its respective group. The order the specimens are sampled determines the number of data points to be removed (i.e. the first to be sampled has the most removed). A complete mathematical description may be found in Brown et al. (In Press).
}
\usage{
byclade(x, remperc, ngroups, groups)
}
\arguments{
  \item{x}{
A n X m matrix of morphometric data with n specimens and m variables
}
  \item{remperc}{
The percentage of data to be removed from the matrix, expressed as a decimal (ex: 30 percent would be entered as 0.3)
}
  \item{ngroups}{
The number of taxonomic groups present in the data matrix
}
  \item{groups}{
A vector of length n specifying taxonomic group membership as integers (ex: c(1,1,2,2,3,3,...) )
}
}
\value{
returns a n X m matrix of morphometric data with missing variables input as 'NA'
}

\references{
Brown, C., Arbour, J. and Jackson, D. 2012. Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology} 61(6):941-954.
}
\author{
J. Arbour and C. Brown
}
\seealso{
\code{\link{missing.data}},\code{\link{obliterator}}
}
\examples{

data(crocs)

## make some fake groups
rows<-nrow(crocs)
groups<-rep(1:4,length.out=rows)
groups

## remove 20% of data points
croc.miss<-byclade(crocs,0.2,4,groups)
croc.miss
}