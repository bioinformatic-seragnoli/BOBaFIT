\name{plot_clust_sample_rep}
\alias{plot_clust_sample_rep}
\docType{data}
\title{
PlotCluster sample clustering report
}
\description{
For each sample, the function produces a dataframe containing the cluster number and the cn of each chromosomal arm
}
\usage{data("plot_clust_sample_rep")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{chr}}{a character vector with the chromosomal arm}
    \item{\code{cluster}}{a character vector with the arm cluster}
    \item{\code{CN}}{a numeric vector with arm Copy Number}
  }
}

\examples{
data(plot_clust_sample_rep)
str(plot_clust_sample_rep)
}
\keyword{datasets}
