\name{align.missing}
\alias{align.missing}
\title{
Procrustes superimposition of landmark datasets with some missing values
}
\description{
This function carries out a generalized procrustes superimposition on all fully complete specimens and produces a consensus confirguration (using "Shapes" procGPA). Each incomplete specimen is then individualyl rotated and aligned with the consensus configuration based on whatever landmarks are available (using "Shapes" procOPA). Landmarks are returned in the same shape space as the original dataset.
}
\usage{
align.missing(X, nlandmarks)
}
\arguments{
  \item{X}{
A n* l X 2 matrix of coordinate data, where n is the number of specimens and l is the number of landmarks. All landmarks from one specimen should be grouped together. 
}
  \item{nlandmarks}{
The number of landmarks per specimen
}
}

\value{
Returns and n*l X 2 matrix of aligned landmarks
}
\references{
Arbour, J. and Brown, C. \emph{In Press}. Incomplete specimens in Geometric Morphometric Analyses. \emph{Methods in Ecology and Evolution}
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
dac.miss<-missing.specimens(dacrya,10,c(1,2,3,4,5,6),16)
dac.miss

## align all specimens
dac.aligned<-align.missing(dac.miss,16)

## compare original and aligned
plot(dacrya,col="blue")
points(dac.aligned,col="red")
}