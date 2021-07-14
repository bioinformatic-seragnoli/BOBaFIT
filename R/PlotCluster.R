#' PlotCluster
#'
#'
#'@description The function clusters chromosomes based on the copy number (CN) and returns a graph where it is possible to observe the different groups and two data frames (report and plot_table). See the vignette for the data frame descriptions.
#'
#' @param segs data.frame with segments of samples. It must be formatted with correct column names (start, end, ID)
#' @param clust_method clustering method. Default is "ward.D2"
#' @param plot_output Whether to plot refitted profiles (logical)
#' @param plot_path Path to save output plots
#'
#' @return Plot with chromosomes clustered
#' @export
#'
#' @importFrom dplyr filter group_by summarise arrange
#' @import NbClust
#' @importFrom ggplot2 ggplot geom_hline geom_point geom_label ggtitle aes
#' @importFrom ggforce geom_mark_ellipse
#' @importFrom stats median weighted.mean
#' @importFrom grDevices png dev.off
#' @importFrom tidyr %>%
#' @importFrom stringr str_sort
#'
#' @examples
#' data(segments)
#' Cluster <- PlotCluster(segs=segments, clust_method= "ward.D2", plot_output=FALSE)
#' Cluster$report
#' Cluster$plot_table

PlotCluster <- function(segs,
                        clust_method = "ward.D2",
                        plot_output= FALSE,
                        plot_path) {
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
          stringsAsFactors = FALSE
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

      if (plot_output == TRUE) {

      png(
        paste0(plot_path, samples[i], "_PlotCluster.png"),
        width = 16,
        height = 4,
        units = "in",
        res = 300
      )

      print(
        ggplot(CLUST_TABLE, aes(
          x = 1:nrow(CLUST_TABLE),
          y = CN,
          colour = cluster
        )) +
          ylim(0,5) +
          geom_mark_ellipse(aes(fill = cluster)) +
          geom_hline(yintercept = 2, alpha = 0.5) +
          geom_point(size = 2) +
          geom_label(aes(label = chr), nudge_y = 0.1) +
          ggtitle(samples[i])
      )
      options(ggplot2. = "viridis")
      dev.off()
      }

    }


    report_clustering <- rbind(report_clustering, samp_report)

    CLUST_TABLE_LIST [[samples[i]]] <- CLUST_TABLE

  }
  OUTPUT <-
    list(report = report_clustering , plot_tables = CLUST_TABLE_LIST)
  OUTPUT
}
