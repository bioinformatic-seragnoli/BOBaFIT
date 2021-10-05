\name{segments}
\alias{segments}
\docType{data}
\title{
Test dataframe
}
\description{
It is an example of BOBaFIT input. The columns name have to be setted with the following names
}
\usage{data("segments")}
\format{
  A data frame with 55778 observations on the following 7 variables.
  \describe{
    \item{\code{chr}}{a numeric vector with segment chromosome}
    \item{\code{start}}{a numeric vector with the segment start }
    \item{\code{end}}{a numeric vector with the segment end}
    \item{\code{arm}}{a character vector chromosomal arm of the segment}
    \item{\code{CN}}{a numeric vector the segment Copy Number}
    \item{\code{ID}}{a character vector with the sample name}
    \item{\code{width}}{a numeric vector the segment width}
  }
}

\examples{
data(segments)
head(segments)
}
\keyword{datasets}
