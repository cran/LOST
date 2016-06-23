\name{complete.specimens}
\alias{complete.specimens}

\title{
Remove incomplete specimens from a landmark dataset
}
\description{
This function takes a dataset containing both complete and incomplete specimens and removes all incomplete specimens. 
}
\usage{
complete.specimens(dataset, nlandmarks)
}
\arguments{
  \item{dataset}{
A n* l X 2 matrix of coordinate data, where n is the number of specimens and l is the number of landmarks. All landmarks from one specimen should be grouped together. 
}
  \item{nlandmarks}{
The number of landmarks per specimen
}
}
\value{
Returns an c * l X 2 matrix of landmark data, where c is the number of complete specimens and l is the number of landmarks.
}
\references{
Arbour, J. and Brown, C. \emph{In Press}. Incomplete specimens in Geometric Morphometric Analyses. \emph{Methods in Ecology and Evolution}
}
\author{
J. Arbour
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{align.missing}},code{\link{MissingGeoMorph}}
}