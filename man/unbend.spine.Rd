\name{unbend.spine}
\alias{unbend.spine}
\title{
Correct for lateral bending in fish geometric morphometric landmarks
}
\description{
Correct for the impact of lateral bending along the spine of a fish in geometric morphometric landmarks. Fits a polynomial function along the length and width of the specimen, determines the perpendicular residuals and arc length along the polynomial and these are used as the new length and width landmarks. Landmarks are first centered and bilaterally aligned using \code{\link{bilat.align}}.

}
\usage{
unbend.spine(coords, land.pairs, deg = 3, restricted = NULL)
}
\arguments{
  \item{coords}{
A matrix of landmark coordinate data. Columns should be coordinates, and rows landmarks.
}
  \item{land.pairs}{
A 2-column matrix giving the bilaterally paired landmarks. One column should be all "left" landmarks and one all "right" landmarks.
}
  \item{deg}{
The degrees for the polynomial function, passed to the function "poly". Typically 2 or 3.
}
  \item{restricted}{
A limited set of landmarks (row numbers for the coords matrix) to use for bilateral alignment. Typically those representing a rigid/fixed structure (e.g., head). Passed to bilat.align.
}
}
\details{
Resulting landmark data is in the same scale as the original landmark configuration. Can be applied over multiple specimens using for-loops or apply functions.
}
\value{
\item{bilat.aligned}{
Provides the bilaterally aligned landmark data as a matrix
}
\item{unbent}{
Provides the bilaterally aligned and unbent landmark data as a matrix
}
}
\references{
Arbour,J.H. In Prep. Get Unbent! R Tools for the removal of arching and bending of fish specimens in geometric morphometric shape analysis
}
\author{
J.H. Arbour
}
\seealso{
\code{\link{bilat.align}}, \code{\link{unbend.tps.poly}}
}
\examples{

\donttest{
data(darters)
library(rgl)
## bilaterally aligned using only head landmarks
lands.unbent<-unbend.spine(darters$coords[,,2],
darters$land.pairs,deg=3, restricted=darters$restricted)$unbent

plot3d(lands.unbent, aspect=FALSE)

}
}
