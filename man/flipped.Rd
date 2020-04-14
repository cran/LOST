\name{flipped}
\alias{flipped}
\title{
Reflected Relabelling
}
\description{
This function carries out reflected relabelling to estimate missing geometric morphometric landmarks using bilateral symmetry following Gunz et al 2009.

A set of 3D landmarks are mirrored and aligned with the original data (using procOPA from package "shapes"). Missing landmarks are interpolated from the mirrored specimen.
}
\usage{
flipped(specimen, land.pairs, show.plot = FALSE, axis = 1)
}
\arguments{
  \item{specimen}{
An l X 3 matrix of coordinate data, where l is the number of landmarks. Some data should be missing and designated with NA.
}
  \item{land.pairs}{
A 2 column matrix, each row should contain row numbers (from matrix specimen) indicating bilateral pairs of landmarks. Unpaired landmarks do not need to be included. See also bilateral symmetry analyses in package "geomorph".
}
  \item{show.plot}{
Optionally plot the specimen using plot3d from rgl. Estimated landmarks are given in red. Defaults to FALSE.
}
  \item{axis}{
Which axis should be mirrored across. Default is x (1).
}
}


\value{
Returns a l X 3 matrix of landmarks.
}
\references{
Gunz P., Mitteroecker P., Neubauer S., Weber G., Bookstein F. 2009. Principles for the virtual reconstruction of hominin crania. Journal of Human Evolution 57:48-62.}
\author{
J. Arbour
}

\seealso{
\code{\link{MissingGeoMorph}}
}
