\name{align.missing}
\alias{align.missing}
\title{
Procrustes superimposition of landmark datasets with some missing values
}
\description{
This function carries out a generalized procrustes superimposition on all fully complete specimens and produces a consensus configuration (using "Shapes" procGPA). Each incomplete specimen is then individually rotated and aligned with the consensus configuration based on any landmarks are available (using "Shapes" procOPA). Data is returned superimposed.}
\usage{
align.missing(X)
}
\arguments{
  \item{X}{
An l X 2 (or 3) X n array of coordinate data, where n is the number of specimens and l is the number of landmarks.
 }
}

\value{
Returns An l X 2 (or 3) X n array of coordinate data
}
\references{
Arbour, J. and Brown, C. 2014. Incomplete specimens in Geometric Morphometric Analyses. \emph{Methods in Ecology and Evolution} 5(1):16-26.
}
\author{
J. Arbour
}

\seealso{
\code{\link{MissingGeoMorph}}
}
\examples{

data(dacrya)

## make some specimens incomplete
dac.miss<-missing.data(dacrya,remsp=0.2,land.vec=c(1,2,3,4,5,6))

## align all specimens
dac.aligned<-align.missing(dac.miss)

}
