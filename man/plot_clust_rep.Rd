\name{plot_clust_rep}
\alias{plot_clust_rep}
\docType{data}
\title{
PlotCluster report
}
\description{
The dataframe reports che clustering outcome and the number of cluster found on each sample
}
\usage{data("plot_clust_rep")}
\format{
  A data frame with 4 observations on the following 3 variables.
  \describe{
    \item{\code{sample}}{a character vector with the sample name}
    \item{\code{clustering}}{a character vector with the cluster outcome (FAIL or SUCCEDED)}
    \item{\code{num_clust}}{a numeric vector with the number of cluster founded}
  }
}

\examples{
data(plot_clust_rep)
head(plot_clust_rep)
}
\keyword{datasets}
