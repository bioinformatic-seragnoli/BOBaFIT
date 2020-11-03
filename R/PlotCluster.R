#' PlotCluster
#'
#'
#'@description The function clusters chromosomes based on the copy number (CN) and returns a graph where it is possible to observe the different groups
#'
#' @param segs data.frame with segments of samples. It must be formatted with correct column names (start, end, ID)
#' @param clust_method clustering method. Default is "ward.D2"
#' @param plot_path Path to save output plots
#'
#' @return Plot with chromosomes clustered
#' @export
#'
#' @importFrom dplyr filter group_by summarise arrange
#' @import NbClust
#' @import ggplot2
#' @importFrom ggforce geom_mark_ellipse
#' @importFrom stats median weighted.mean
#' @importFrom grDevices png dev.off
#' @importFrom tidyr %>%
#'
#' @examples
#' \dontrun{
#' PlotCluster(segments, clust_method=ward.D2, plot_path)
#' }
PlotCluster <- function(segs, clust_method, plot_path) {
  report_clustering <- data.frame(sample = character(),
                                  clustering = character(),
                                  num_clust = numeric())
  ID <- arm <- CN <- width <- str_sort <- chr <- cluster <- NULL

  CLUST_TABLE_LIST <- list()


  samples <- segs$ID %>% unique()

  for (i in 1:length(samples)) {
    cat("sample n: ", i, " - ", samples[i], "\n")

    segments <- segs %>%  filter(ID == samples[i])

    CN_CHR <- segments %>%
      group_by(arm) %>%
      summarise(weighted_mean_CN = weighted.mean(CN, w = width),
                .groups = 'drop')

    CN_CHR_values <- CN_CHR$weighted_mean_CN


    TRY <- try({
      ClustRes <-
        NbClust(
          CN_CHR_values,
          distance = "euclidean",
          method = clust_method,
          index = "all",
          min.nc = 2,
          max.nc = 6
        )
      CLUST_TABLE <-
        data.frame(
          chr = CN_CHR$arm,
          cluster = ClustRes$Best.partition,
          CN = CN_CHR_values,
          stringsAsFactors = F
        )

    })


    if (class(TRY) == "try-error") {
      print("Clustering failed")
      samp_report <- data.frame(sample = samples[i],
                                clustering = "SUCCEDED",
                                num_clust = NA)

    } else {
      print("Clustering succeded")

      samp_report <- data.frame(
        sample = samples[i],
        clustering = "SUCCEDED",
        num_clust = max(ClustRes$Best.partition)
      )


      CLUST_TABLE$chr <-
        CLUST_TABLE$chr %>% factor(levels = str_sort(CLUST_TABLE$chr, numeric = TRUE),
                                   ordered = TRUE)

      CLUST_TABLE$cluster <- paste0("cluster", CLUST_TABLE$cluster)

      CLUST_TABLE <- CLUST_TABLE %>% arrange(chr)

      png(
        paste0(plot_path, samples[i], "_PlotCluster.png"),
        width = 16,
        height = 4,
        units = "in",
        res = 300
      )

      print(
        ggplot(CLUST_TABLE, aes(
          x = 1:40,
          y = CN,
          colour = cluster
        )) +
          geom_mark_ellipse(aes(fill = cluster)) +
          geom_hline(yintercept = 2, alpha = 0.5) +
          geom_point(size = 2) +
          geom_label(aes(label = chr), nudge_y = 0.1) +
          ggtitle(samples[i])
      )
      options(ggplot2. = "viridis")
      dev.off()
    }


    report_clustering <- rbind(report_clustering, samp_report)

    CLUST_TABLE_LIST [[samples[i]]] <- CLUST_TABLE

  }
  OUTPUT <-
    list(report = report_clustering , plot_tables = CLUST_TABLE_LIST)

}
