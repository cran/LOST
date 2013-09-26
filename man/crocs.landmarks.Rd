\name{crocs.landmarks}
\alias{crocs.landmarks}
\docType{data}
\title{
Coordinate data for a crocodilian reference skull
}
\description{
Landmark data for the measurements points on a reference crocodilian skull, for use with the \code{\link{obliterator}} function
}
\usage{data(crocs.landmarks)}
\format{
  A 6 X m dataframe in which each column gives the start and end points for each cranial measurement in the crocs dataset, from a single reference specimen. 3D Coordinates are listed as x1, x2, y1, y2, z1, z2 in each column.
}
\source{
Brown, C., Arbour, J. and Jackson, D. 2012. Testing of the Effect of Missing Data Estimation and Distribution in Morphometric Multivariate Data Analyses. \emph{Systematic Biology} 61(6):941-954.
}
\seealso{
\code{\link{obliterator}}, \code{\link{byclade}},\code{\link{missing.data}},\code{\link{crocs}}
}