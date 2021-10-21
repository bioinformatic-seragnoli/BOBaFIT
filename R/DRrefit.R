#' DRrefit
#'
#' @description This function refits the diploid region of input copy number profiles (segments - BED file)
#'
#' @param segments_chort data.frame formatted with correct column names
#' @param chrlist list of normal chromosome arms (pathology-specific)
#' @param maxCN threshold of max copy number to consider. By default is 6
#' @param clust_method clustering method. By default is "ward.D2"

#' @return Return two data frames, one is the DRrefit-corrected segments and the other is the samples report. See the vignette for data frame descriptions.
#' @export
#'
#'
#' @importFrom dplyr filter group_by summarise arrange n desc
#' @import NbClust
#' @importFrom stats median weighted.mean
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tidyr %>%
#'
#' @examples
#' segments <- data("TCGA_BRCA_CN_segments")
#' chr_list <- c("10q","11p","12p","19q","1p","21q","2q","3p","4p","4q","6p","6q","7p" )
#' 
#' results <- DRrefit(segments,chrlist = chr_list)
#' results$report
#' results$segments_corrected
TCGA_BRCA_CN_segments <- segments
DRrefit <- function(segments_chort,
                    chrlist,
                    maxCN = 6,
                    clust_method = "ward.D2") {

  ward.D2 <- ID <- chrarm <- CN <- width <- chr <- cluster <- CN_corrected <- number <- NULL
  OUTPUT <- list()

  segments_chort$CN[segments_chort$CN == Inf] <- maxCN

  samples <- segments_chort$ID %>% unique()

  segments_chort_corrected <- data.frame()
  report_clustering <- data.frame(sample=character(),
                                  clustering=character(),
                                  ref_clust_chr=character(),
                                  num_clust=numeric(),
                                  correction_factor=numeric(),
                                  correction_class= character())

  for (i in seq_along(samples)){

    cat("sample n: ",i," - ",samples[i],"\n")

    segments <- segments_chort %>%  filter(ID==samples[i])

    CN_CHR <- segments %>%
      group_by(chrarm) %>%
      summarise(weighted_mean_CN = weighted.mean(CN, w = width), .groups = 'drop')


    CN_CHR_values <- CN_CHR$weighted_mean_CN


    TRY <- try({


      ClustRes <- NbClust(CN_CHR_values, distance = "euclidean", method = clust_method, index = "all", min.nc = 2, max.nc = 6)


      CLUST_TABLE <- data.frame(chr=CN_CHR$chrarm, cluster=ClustRes$Best.partition, stringsAsFactors = FALSE)


      CLUST_TABLE_in_list_sorted <- CLUST_TABLE %>%
        filter(chr %in% chrlist ) %>%
        group_by(cluster) %>%
        summarise(number = n()) %>%
        arrange(desc(number))


      ref_cluster <- CLUST_TABLE_in_list_sorted$cluster[1]



      ref_cluster <- CLUST_TABLE_in_list_sorted$cluster[1]


      ref_cluster_chrs <- CLUST_TABLE %>% filter(cluster==ref_cluster)
      chrlist_cluster <- ref_cluster_chrs$chr

    })

    if (is(TRY, "try-error")){
      message("Clustering failed")

      new_chrlist <-  chrlist

      samp_report <-data.frame(
        sample=samples[i],
        clustering="FAIL",
        ref_clust_chr=paste(chrlist, collapse = ", "),
        num_clust=NA)




    } else {
      message("Clustering succeded")

      new_chrlist <-  chrlist_cluster

      samp_report <-data.frame(
        sample=samples[i],
        clustering="SUCCEDED",
        ref_clust_chr=paste(chrlist_cluster, collapse = ", "),
        num_clust=max(ClustRes$Best.partition))


    }


    segments_REF <- segments %>% filter(chrarm %in% new_chrlist)

    segments_REF <- segments %>% filter(chrarm %in% new_chrlist)

    CN_CHR_REF <- segments_REF %>%
      group_by(chrarm) %>%
      summarise(weighted_mean_CN = weighted.mean(CN, w = width), .groups = 'drop')

    real_diploid_region <-  median(CN_CHR_REF$weighted_mean_CN)

    correction_factor <- 2 - real_diploid_region

    samp_report$correction_factor <- correction_factor

    samp_report$correction_class <- ifelse(correction_factor > 0.5, "REFITTED",
                                           ifelse(correction_factor <= 0.1, "NO CHANGES", "RECALIBRATED"))

    report_clustering <- rbind(report_clustering, samp_report)

    segments$CN_corrected <- segments$CN + correction_factor

    segments$CN_corrected <- ifelse(segments$CN_corrected < 0, 0.001, segments$CN_corrected)


    segments_chort_corrected <- rbind(segments_chort_corrected, segments)
  }

  segments_chort_corrected

  OUTPUT[["segments_corrected"]] <- segments_chort_corrected
  OUTPUT[["report"]] <- report_clustering

  OUTPUT
}
