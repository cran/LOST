\name{unbend.tps.poly}
\alias{unbend.tps.poly}

\title{
TPS-style unbend specimens
}
\description{
Remove dorsoventral arching effect from fish specimen landmark data. Function similar to "unbend specimens" utility in the TPS software suite. Fits a polynomial function along the length and height of the specimen, determines the perpendicular residuals and arc length along the polynomial, and these are used as the new length and width landmarks. 
}
\usage{
unbend.tps.poly(coords, reference, axes = NULL, deg = 3)
}
\arguments{
  \item{coords}{
A matrix of landmark coordinate data. Columns should be coordinates, and rows landmarks.
}
  \item{reference}{
The rows of the matrix over which the polynomial function will be fit. Should represent the spine or other proxy for the long axis of the body.}
  \item{axes}{
A vector with 2 values representing the "lateral" view of the fish. The first entry should be the "long" (anterior-posterior) axis and the second should be the vertical (dorso-central) axis.
}
  \item{deg}{
The degrees for the polynomial function, passed to "poly". Typically 2 or 3 (default = 3).
}
}
\details{
It is advisable to remove lateral bending with \code{\link{unbend.spine}} prior to using this function. Otherwise data should be at least bilaterally aligned to a plane (and see\code{\link{bilat.align}}) Resulting landmark data is in the same scale as the original landmark configuration. Can be applied over multiple specimens using for-loops or apply functions.
}
\value{
Returns a matrix of landmark data with the effect of dorso-ventral arching removed. 
}
\references{
Arbour,J.H. In Prep. Get Unbent! R Tools for the removal of arching and bending of fish specimens in geometric morphometric shape analysis
}
\author{
J.H. Arbour
}

\seealso{
\code{\link{bilat.align}}, \code{\link{unbend.spine}}
}
\examples{
\donttest{
library(rgl)
data(darters)
## bilaterally aligned using only head landmarks
lands.unbent<-unbend.spine(darters$coords[,,3],
darters$land.pairs,deg=3, restricted=darters$restricted)$unbent

plot(lands.unbent[,c(1,3)],asp=1)

lands.unbent<-unbend.tps.poly(lands.unbent,darters$reference,axes=c(1,3))
plot(lands.unbent[,c(1,2)],asp=1)

plot3d(lands.unbent, aspect=FALSE)
}
}