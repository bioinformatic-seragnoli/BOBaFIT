#' PlotChrCluster
#'
#'
#'@description The function clusters chromosomes based on the copy number (CN) and returns a graph where it is possible to observe the different groups and two data frames (report and plot_table). See the vignette for the data frame descriptions.
#'
#' @param segs data.frame with segments of samples. It must be formatted with correct column names (start, end, ID)
#' @param clust_method clustering method. Default is "ward.D2"
#' @param plot_output Whether to plot refitted profiles (logical)
#' @param plot_viewer Logical parameter. When it is TRUE, the function print the output plot in the R viewer.By default is TRUE.
#' @param plot_save  Logical parameter. When it is TRUE, the function save the plot in the chosen path and format. By default is TRUE.
#' @param plot_format File format for the output plots (accepts "png", "jpg", "pdf", "tiff"). By default is "png" 
#' @param plot_path Path to save output plots.
#' @param verbose print information about the processes of the function. By default is FALSE
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
#' @importFrom stringr str_sort
#' @importFrom methods is
#' @importFrom utils capture.output
#'
#' @examples
#' data(TCGA_BRCA_CN_segments)
#' Cluster <- PlotChrCluster(segs=TCGA_BRCA_CN_segments, 
#'                          clust_method= "ward.D2", 
#'                          plot_output=FALSE)

PlotChrCluster <- function(segs,
                           clust_method = "ward.D2",
                           plot_output= TRUE,
                           plot_viewer= TRUE,
                           plot_save = FALSE,
                           plot_format = "png",
                           plot_path,
                           verbose= FALSE) {
  report_clustering <- data.frame(sample = character(),
                                  clustering = character(),
                                  num_clust = numeric())
  ID <- chrarm <- CN <- width <- str_sort <- chr <- cluster <- NULL
  
  CLUST_TABLE_LIST <- list()
  
  
  samples <- segs$ID %>% unique()
  if (verbose) {
    for (i in seq_along(samples)) {
      cat("sample n: ", i, " - ", samples[i], "\n")
      
      segments <- segs %>%  filter(ID == samples[i])
      
      CN_CHR <- segments %>%
        group_by(chrarm) %>%
        summarise(weighted_mean_CN = weighted.mean(CN, w = width),
                  .groups = 'drop')
      
      CN_CHR_values <- CN_CHR$weighted_mean_CN
      
      
      TRY <- try({
        capture.output(ClustRes <- NbClust(CN_CHR_values, distance = "euclidean", method = clust_method, index = "all", min.nc = 2, max.nc = 6),  type = c("output", "message"))
        
        CLUST_TABLE <-
          data.frame(
            chr = CN_CHR$chrarm,
            cluster = ClustRes$Best.partition,
            CN = CN_CHR_values,
            stringsAsFactors = FALSE
            
          )
      }, silent = TRUE)
      
      
      if (is(TRY, "try-error")) {
        message("Clustering failed")
        samp_report <- data.frame(sample = samples[i],
                                  clustering = "FAIL",
                                  num_clust = NA)
        
      } else {
        message("Clustering succeded")
        
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
        
        if (plot_output) {
          
          gg <- ggplot(CLUST_TABLE, aes(
            x = seq_len(nrow(CLUST_TABLE)),
            y = CN, colour=cluster)) +
            ylim(-0.5, max(CLUST_TABLE$CN)) +
            geom_mark_ellipse( aes(fill = cluster)) +
            geom_hline(yintercept = 2, alpha = 0.5) +
            geom_point(size = 2) +
            geom_label(aes(label = chr), nudge_y = 0.1) +
            ggtitle(samples[i])+
            xlab("Chromosomal arms")+ 
            ylab("Copy Number")
          
          options(ggplot2. = "viridis")
          
          if(plot_save){
            
            if ( plot_format=="png" ) { png( paste0(plot_path, samples[i],"_PlotChrCluster.png"), width = 16, height = 4, units = "in", res = 300 ) 
            } else if(plot_format=="tiff") { tiff( paste0(plot_path, samples[i],"_PlotChrCluster.tif"), width = 16, height = 4, units = "in", res = 300 )
            } else if(plot_format=="pdf") { pdf( file = paste0(plot_path, samples[i],"_PlotChrCluster.pdf"), width = 16, height = 4 )
            } else if(plot_format == "jpg") { jpeg( paste0(plot_path, samples[i],"_PlotChrCluster.jpg"), width = 16, height = 4, units = "in", res = 300 ) 
            } else { message("wrong plot format") 
              break}
            
            print(gg)
            
            dev.off()
          }
          
          if(plot_viewer){
            print(gg)
          }
          
        }
        
        
        report_clustering <- rbind(report_clustering, samp_report)
        
        CLUST_TABLE_LIST [[samples[i]]] <- CLUST_TABLE
        
      }
    }} else {
      
      for (i in seq_along(samples)) {
        
        segments <- segs %>%  filter(ID == samples[i])
        
        CN_CHR <- segments %>%
          group_by(chrarm) %>%
          summarise(weighted_mean_CN = weighted.mean(CN, w = width),
                    .groups = 'drop')
        
        CN_CHR_values <- CN_CHR$weighted_mean_CN
        
        
        TRY <- try({
          capture.output(ClustRes <- NbClust(CN_CHR_values, distance = "euclidean", method = clust_method, index = "all", min.nc = 2, max.nc = 6),  type = c("output", "message"))
          
          CLUST_TABLE <-
            data.frame(
              chr = CN_CHR$chrarm,
              cluster = ClustRes$Best.partition,
              CN = CN_CHR_values,
              stringsAsFactors = FALSE
              
            )
        }, silent = TRUE)
        
        
        if (is(TRY, "try-error")) {
          samp_report <- data.frame(sample = samples[i],
                                    clustering = "FAIL",
                                    num_clust = NA)
          
        } else {
          
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
          
          if (plot_output) {
            
            gg <- ggplot(CLUST_TABLE, aes(
              x = seq_len(nrow(CLUST_TABLE)),
              y = CN, colour=cluster)) +
              ylim(-0.5, max(CLUST_TABLE$CN)) +
              geom_mark_ellipse( aes(fill = cluster)) +
              geom_hline(yintercept = 2, alpha = 0.5) +
              geom_point(size = 2) +
              geom_label(aes(label = chr), nudge_y = 0.1) +
              ggtitle(samples[i])+
              xlab("Chromosomal arms")+ 
              ylab("Copy Number")
            
            
            options(ggplot2. = "viridis")
            
            if(plot_save){
              
              if ( plot_format=="png" ) { png( paste0(plot_path, samples[i],"_PlotChrCluster.png"), width = 16, height = 4, units = "in", res = 300 ) 
              } else if(plot_format=="tiff") { tiff( paste0(plot_path, samples[i],"_PlotChrCluster.tif"), width = 16, height = 4, units = "in", res = 300 )
              } else if(plot_format=="pdf") { pdf( file = paste0(plot_path, samples[i],"_PlotChrCluster.pdf"), width = 16, height = 4 )
              } else if(plot_format == "jpg") { jpeg( paste0(plot_path, samples[i],"_PlotChrCluster.jpg"), width = 16, height = 4, units = "in", res = 300 ) 
              } else { message("wrong plot format") 
                break}
              
              print(gg)
              
              dev.off()
            }
            
            if(plot_viewer){
              print(gg)
            }
            
          }
          
          
          
          
          report_clustering <- rbind(report_clustering, samp_report)
          
          CLUST_TABLE_LIST [[samples[i]]] <- CLUST_TABLE
          
        }} }
  
  OUTPUT <-
    list(report = report_clustering , plot_tables = CLUST_TABLE_LIST)
  OUTPUT
}
