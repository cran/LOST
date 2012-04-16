\name{probability.generator}
\alias{probability.generator}
\title{
Internal function for obliterator
}
\description{
Internal function
}
\usage{
probability.generator(newsp, distances, removes, expo = 1, sa)
}
\arguments{
  \item{newsp}{
specimen to have data removed from
}
  \item{distances}{
reference landmarks
}
  \item{removes}{
number of measurements to remove
}
  \item{expo}{
an optional exponent for increasing or decreasing the intensity of anatomical bias
}
  \item{sa}{
a vector tracking the removed measurements
}
}

\seealso{
\code{\link{obliterator}}
}
