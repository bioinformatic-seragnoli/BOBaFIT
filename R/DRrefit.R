#' DRrefit
#'
#' @description This function refits the diploid region of input copy number profiles (BED file)
#'
#' @param segments_chort data.frame formatted with correct column names
#' @param chrlist list of normal chromosome arms (pathology-specific)
#' @param maxCN threshold of max copy number to consider. By default is 6
#' @param clust_method clustering method. By default is "ward.D2"
#' @param plot_output Whether to plot refitted profiles (logical)
#' @param plot_path Path to save output plots
#'
#' @return A plot that shows old and refitted copy number profiles  and two databases (segment corrected and report)
#' @export
#'
#'
#' @importFrom dplyr filter group_by summarise arrange n desc
#' @import NbClust
#' @importFrom ggplot2 ggplot ggtitle ylim facet_grid theme_bw theme scale_x_continuous geom_hline  xlab scale_color_manual unit element_rect
#' @importFrom ggbio geom_segment
#' @importFrom stats median weighted.mean
#' @importFrom grDevices png dev.off
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tidyr %>%
#'
#' @examples
#'\donttest{
#' data(segments)
#' chrlist <- c("10q","11p","12p",19q","1p","21q","2q","3p","4p","4q","6p","6q","7p" )
#' DRrefit(segments,chrlist = chrlist , maxCN=6, clust_method= "ward.D2", plot_output=FALSE)
#' }

DRrefit <- function(segments_chort,
                    chrlist,
                    maxCN = 6,
                    clust_method = "ward.D2",
                    plot_output = FALSE,
                    plot_path) {

  ward.D2 <- ID <- arm <- CN <- width <- chr <- cluster <- CN_corrected <- number <- NULL
  OUTPUT <- list()

  segments_chort$CN[segments_chort$CN == Inf] <- maxCN

  segments_chort$width <- segments_chort$end - segments_chort$start

  samples <- segments_chort$ID %>% unique()

  segments_chort_corrected <- data.frame()
  report_clustering <- data.frame(sample=character(),
                                  clustering=character(),
                                  ref_clust_chr=character(),
                                  num_clust=numeric(),
                                  correction_factor=numeric())

  for (i in 1:length(samples)){

    cat("sample n: ",i," - ",samples[i],"\n")

    segments <- segments_chort %>%  filter(ID==samples[i])

    CN_CHR <- segments %>%
      group_by(arm) %>%
      summarise(weighted_mean_CN = weighted.mean(CN, w = width), .groups = 'drop')


    CN_CHR_values <- CN_CHR$weighted_mean_CN


    TRY <- try({


      ClustRes <- NbClust(CN_CHR_values, distance = "euclidean", method = clust_method, index = "all", min.nc = 2, max.nc = 6)


      CLUST_TABLE <- data.frame(chr=CN_CHR$arm, cluster=ClustRes$Best.partition, stringsAsFactors = FALSE)


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

    if (class(TRY)=="try-error"){
      print("Clustering failed")

      new_chrlist <-  chrlist

      samp_report <-data.frame(
        sample=samples[i],
        clustering="FAIL",
        ref_clust_chr=paste(chrlist, collapse = ", "),
        num_clust=NA)




    } else {
      print("Clustering succeded")

      new_chrlist <-  chrlist_cluster

      samp_report <-data.frame(
        sample=samples[i],
        clustering="SUCCEDED",
        ref_clust_chr=paste(chrlist_cluster, collapse = ", "),
        num_clust=max(ClustRes$Best.partition))


    }


    segments_REF <- segments %>% filter(arm %in% new_chrlist)

    segments_REF <- segments %>% filter(arm %in% new_chrlist)

    CN_CHR_REF <- segments_REF %>%
      group_by(arm) %>%
      summarise(weighted_mean_CN = weighted.mean(CN, w = width), .groups = 'drop')

    real_diploid_region <-  median(CN_CHR_REF$weighted_mean_CN)

    correction_factor <- 2 - real_diploid_region

    samp_report$correction_factor <- correction_factor

    report_clustering <- rbind(report_clustering, samp_report)

    segments$CN_corrected <- segments$CN + correction_factor

    segments$CN_corrected <- ifelse(segments$CN_corrected < 0, 0.001, segments$CN_corrected)

    if (plot_output == TRUE) {

      Granges_segments <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)

      allchr <- segments$arm %>% unique()
      idx <- allchr %in% new_chrlist

      altered_chr <- allchr[!idx]

      png( paste0(plot_path, samples[i],"_",clust_method,".png"), width = 16, height = 4, units = "in", res = 300 )

      print(
        ggplot(Granges_segments) +
          ggtitle(paste0("Sample name: ",samples[i]),
                  subtitle = paste0("Correction factor= ", correction_factor %>% round(3)," / not diploid chrs: ", paste(altered_chr, collapse = ", ")) ) +
          geom_segment(stat="identity", aes(y=CN, colour="old_CN"), size = 1.2, alpha=0.7) +
          geom_segment(stat="identity", aes(y=CN_corrected, colour="new_CN_corrected"), size = 1.2, alpha=0.7) +
          ylim(0,5) +
          facet_grid( ~seqnames, scales = "free", space = "free", margins = FALSE) +
          theme_bw() +
          theme(panel.spacing=unit(.05, "lines"),
                panel.border = element_rect(color = "black", fill = NA, size = 0.2),
                legend.position="bottom") +
          scale_x_continuous(labels = function(x) format(x/1000000, scientific = FALSE), expand = c(0, 0), limits = c(0, NA))+
          geom_hline(yintercept = 2, linetype=3) +
          xlab("Mbp (genomic position)") +
          if(abs(correction_factor) > 0.1 ){
            scale_color_manual(values=c("green4", "red3"))
          } else {
            scale_color_manual(values=c("orange", "red3"))
          }
      )

      options(ggplot2. ="viridis")

      dev.off()

    }


    segments_chort_corrected <- rbind(segments_chort_corrected, segments)
  }

  segments_chort_corrected

  OUTPUT[["segments_corrected"]] <- segments_chort_corrected
  OUTPUT[["report"]] <- report_clustering

  OUTPUT
}
