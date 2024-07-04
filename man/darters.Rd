\name{darters}
\alias{darters}
\docType{data}
\title{
Darter landmarks
}
\description{
A 3D landmark dataset from 30 species of darter fishes (Etheostomatinae; Percidae)
}
\usage{data("darters")}
\format{
  The format is:
List of 6
 $ coords    : num [1:220, 1:3, 1:30] -1.458 -0.489 -0.037 1.705 0.959 ...
  ..- attr(*, "dimnames")=List of 3
  .. ..$ : NULL
  .. ..$ : NULL
  .. ..$ : chr [1:30] "Etheostoma_caeruleum_mtsu5_58mmsl.stl" "Ammocrypta_beanii_ummz242736_43mm.stl" "Ammocrypta_clara_ummz148570_42.23mm.stl" "Crystallaria_asprella_Ummz211889_60mmSL.stl" ...
 $ land.pairs:'data.frame':	101 obs. of  2 variables:
  ..$ left : int [1:101] 1 3 5 7 9 11 13 15 17 19 ...
  ..$ right: int [1:101] 2 4 6 8 10 12 14 16 18 20 ...
 $ sliders   :'data.frame':	32 obs. of  3 variables:
  ..$ start: int [1:32] 22 23 24 25 26 27 28 29 31 32 ...
  ..$ slide: int [1:32] 23 24 25 26 27 28 29 30 32 33 ...
  ..$ end  : int [1:32] 24 25 26 27 28 29 30 31 33 34 ...
 $ surface   :'data.frame':	144 obs. of  1 variable:
  ..$ surface: int [1:144] 60 61 62 63 64 65 66 68 69 70 ...
 $ restricted: int [1:58] 1 2 3 4 5 6 7 8 9 10 ...
 $ reference : num [1:11] 22 99 180 15 16 63 176 81 178 11 ...
}
\details{
Includes landmark coordinates (coords), a matrix indicating bilaterally paired landmarks (land.pairs), curve sliders (sliders), surface sliders (surface), rows of head landmarks (restricted) and landmarks approximating the spine/long axis (reference). 
}
\source{
Arbour,J.H. In Prep. Get Unbent! R Tools for the removal of arching and bending of fish specimens in geometric morphometric shape analysis
}
\references{
Arbour,J.H. In Prep. Get Unbent! R Tools for the removal of arching and bending of fish specimens in geometric morphometric shape analysis
}
\examples{
data(darters)
library(rgl)
plot3d(darters$coords[,,1], aspect=FALSE)
}
\seealso{
\code{\link{unbend.spine}},\code{\link{bilat.align}},\code{\link{unbend.tps.poly}}
}
\keyword{datasets}
