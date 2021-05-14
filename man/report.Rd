\name{report}
\alias{report}
\docType{data}
\title{
DRrefit cluster report
}
\description{
The dataframe contains all the information about the clustering as the outcome,the number of clusters found in that sample, the chromosome list and the CR used for the adjstment of the diploid region
}
\usage{data("report")}
\format{
  A data frame with 100 observations on the following 5 variables.
  \describe{
    \item{\code{sample}}{a character vector with sample name}
    \item{\code{clustering}}{a character vector with the clustering outcome (FAIL or SUCCEDED)}
    \item{\code{ref_clust_chr}}{a character vector with the chromosme list used}
    \item{\code{num_clust}}{a numeric vector with the number of cluster founded in the sample}
    \item{\code{correction_factor}}{a numeric vector with the correction factor used to move the dipolid region}
  }
}

\examples{
data(report)
str(report)
}
\keyword{datasets}
