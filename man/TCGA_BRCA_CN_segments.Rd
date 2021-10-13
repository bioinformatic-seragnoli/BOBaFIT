\name{TCGA_BRCA_CN_segments.RData}
\alias{TCGA_BRCA_CN_segments.RData}
\docType{data}
\title{
Test dataframe
}
\description{
Dataset downloaded from TCGA databases, it is composed by Breast Cancer (BRCA) samples' segments obtained by SNParray analysis (Affimetrix SNP 6.0 platform and DNAcopy tool). We just added the arm column and in the good practice vignette are indicated the steps to prepare the dataset for BOBaFIT's analysis. 
}
\usage{data("segments")}
\format{
  A data frame with 55778 observations on the following 7 variables.
\describe{
 \item{\code{arm}}{a character vector with the chromosomal arm of segment}
 \item{\code{chr}}{a numeric vector with the chromosome of the segment}
 \item{\code{start}}{a numeric vector with the segment start}
 \item{\code{end}}{a numeric vector with the segment end}
 \item{\code{width}}{a numeric vector with the segment width (Mb)}
 \item{\code{strand}}{a factor vector with the segment strand sense (levels: +,-,*)}
 \item{\code{GDCAliquot}}{a character vector with the sample name}
 \item{\code{Num_Probes}}{a numeric vector with the number of probes}
 \item{\code{Segment_mean}}{a numeirc vector with mthe mean signal of the probes}
 \item{\code{Sample}}{a charcter vector with the TCGA-BRCA barcode}
  }
}

\examples{
data(segments)
head(segments)
}
\keyword{datasets}
