\name{missing.specimens}
\alias{missing.specimens}
\title{
Simulate incomplete specimens
}
\description{
Randomly selects a pre-selected number of specimens from a landmark dataset and removes some of its landmarks.
}
\usage{
missing.specimens(dataset, nspremove, nldremove, nlandmarks)
}
\arguments{
  \item{dataset}{
A n* l X 2 matrix of coordinate data, where n is the number of specimens and l is the number of landmarks. All landmarks from one specimen should be grouped together. 
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
Returns an n * l X 2 matrix with some complete and some incomplete specimens.
}
\references{
Arbour, J. and Brown, C. \emph{In Press}. Incomplete specimens in Geometric Morphometric Analyses. \emph{Methods in Ecology and Evolution}
}
\author{
J. Arbour
}

\seealso{
\code{\link{align.missing}}, \{code{\link{MissingGeoMorph}}
}
\examples{
data(dacrya)

#### remove 1 to 6 landmarks from 10 specimens
dac.miss<-missing.specimens(dacrya,10,c(1,2,3,4,5,6),16)
dac.miss
}
