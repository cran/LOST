\name{missing.specimens}
\alias{missing.specimens}
\title{
Simulate incomplete specimens
}
\description{
Randomly selects a pre-determined number of specimens from a landmark dataset (2D or 3D) and removes some of their landmarks.
}
\usage{
missing.specimens(dataset, nspremove, nldremove, nlandmarks)
}
\arguments{
  \item{dataset}{
A n*l X 2 (or 3) matrix of coordinate data, where n is the number of specimens and l is the number of landmarks. All landmarks from one specimen should be grouped together. 
}
  \item{nspremove}{
The number of specimens which should have landmarks removed. 
}
  \item{nldremove}{
The number of landmarks to remove per specimen. This may be a single value or a vector of values, none of which can be >nlandmarks. If a vector is given, for each specimen selected, the function will randomly select a value from the vector and remove that many landmarks.
}
  \item{nlandmarks}{
The number of landmarks per specimen
}
}
\value{
Returns an n * l X 2 (or 3) matrix with some complete and some incomplete specimens.
}
\references{
Arbour, J. and Brown, C. 2014. Incomplete specimens in Geometric Morphometric Analyses. \emph{Methods in Ecology and Evolution} 5(1):16-26.
}
\author{
J. Arbour
}

\seealso{
\code{\link{align.missing}}, \code{\link{MissingGeoMorph}}
}

