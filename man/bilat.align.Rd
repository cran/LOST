\name{bilat.align}
\alias{bilat.align}
\title{
Align a bilterally symmetric landmark configuration with a plane
}
\description{
Aligns a bilaterally symmetric landmark dataset to a specific plane by minimized the sum of squared distances of one coordinate (x, y or z). Useful for averaging bilateral landmarks or in preparation for correcting for artifacts like bending.

}
\usage{
bilat.align(coords, land.pairs, average = TRUE, restricted = NULL)
}
\arguments{
  \item{coords}{
Either a matrix or array of landmark data with columns representing the x, y, z coordinates and rows representing landmarks. See details for how this is applied for a single vs. multiple specimens.
}
  \item{land.pairs}{
A 2 column matrix indicating bilaterally paired landmarks. All "left" landmarks should be in the same column (and likewise for "right landmarks")
}
  \item{average}{
An optional term indicating that bilaterally paired landmarks should be mirrored and averaged, leaving only one "side" and the midline landmarks.
}
  \item{restricted}{
A set of row numbers indicating which landmarks should be considered by "optim" when selecting the optimal rotation. Typically landmarks representing a rigid structure if some landmarks represent articulated/moveable features.
}
}
\details{
If a matrix for a single specimen's landmarks is provided this is aligned to a plane, if an array of multiple specimens is provided, these should be previously aligned with Procrustes superimposition, and the entire configuration is optimized with a single rotation applied to all specimens. SS are minimized across the third axis (coords[,3] or coords[,3,]).
}
\value{
A matrix or array giving the rotated landmark configuration
}
\references{
Arbour,J.H. In Prep. Get Unbent! R Tools for the removal of arching and bending of fish specimens in geometric morphometric shape analysis
}
\author{
J.H. Arbour
}

\seealso{
\code{\link{unbend.spine}}, \code{\link{unbend.tps.poly}}
}
\examples{
library(rgl)
data(darters)
## align darter configuration by head landmarks (restricted)
aligned<-bilat.align(darters$coords[,,1],
darters$land.pairs,average=FALSE,darters$restricted)

plot3d(aligned, aspect=FALSE)

}