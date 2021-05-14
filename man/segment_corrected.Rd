\name{segment_corrected}
\alias{segment_corrected}
\docType{data}
\title{
Segment corrected by DRrefit
}
\description{
It contains the segment the old and the new CN corrected by the correction factor (CR) of all samples analyzed
}
\usage{data("segment_corrected")}
\format{
  A data frame with 55778 observations on the following 8 variables.
  \describe{
    \item{\code{chr}}{a numeric vector with segment chromosome}
    \item{\code{start}}{a numeric vector with the segment start}
    \item{\code{end}}{a numeric vector with the segment end}
    \item{\code{arm}}{a character vector chromosomal arm of the segment}
    \item{\code{CN}}{a numeric vector the old segment Copy Number}
    \item{\code{ID}}{a character vector with the sample name}
    \item{\code{width}}{a numeric vector the segment width}
    \item{\code{CN_corrected}}{a numeric vector the new segment Copy Number}
  }
}

\examples{
data(segment_corrected)
head(segment_corrected)
}
\keyword{datasets}
