#' DRrefit_plot
#'
#' @param segments_refitted DRrefit output dataframe. 
#' @param DRrefit_report  DRrefit output dataframe.
#' @param plot_viewer Logical parameter. When it is TRUE, the function print the output plot in the R viewer.By default is FALSE.
#' @param plot_save  Logical parameter. When it is TRUE, the function save the plot in the chosen path and format. By default is TRUE.
#' @param plot_format File format for the output plots. By default is "png" (accepts "png", "jpg", "pdf", "tiff").
#' @param plot_path Path to save output plots.
#'
#' @return Return the sample copy number profile before and after DRrefit recalibration. The function can output the figure in the R viewer on save it in a specific path.
#' @export
#' 
#' @importFrom grDevices png tiff pdf jpeg dev.off
#' @importFrom ggplot2 ggplot ggtitle ylim facet_grid theme_bw theme scale_x_continuous geom_hline  xlab scale_color_manual unit element_rect
#' @importFrom ggbio geom_segment
#' @importFrom GenomicRanges makeGRangesFromDataFrame 
#' @importFrom tidyr %>%
#' @importFrom dplyr filter 
#'
#' @examples
#' 
#' segments <- data(TCGA_BRCA_CN_segments)
#' chr_list <- c("10q","11p","12p","19q","1p","21q","2q","3p","4p","4q","6p","6q","7p" )
#' 
#' results <- DRrefit(segments,chrlist = chr_list)
#' DRrefit_report <- results$report
#' segments_refitted <- results$segments_corrected
#' 
#' DRrefit_plot(segments_refitted,DRrefit_report, plot_viewer= TRUE, plot_save = FALSE)
#' 
#' 

DRrefit_plot <- function(segments_refitted,
                         DRrefit_report,
                         plot_viewer=FALSE,
                         plot_save = TRUE,
                         plot_format = "png",
                         plot_path
                         ) {
  
  samples <- unique(segments_refitted$ID)
  
  for (i in seq_along(samples)) {
    
    segments <- segments_refitted %>%  filter(ID==samples[i])
    
    report_samp <- report %>%  filter(sample== samples[i])
    
    Granges_segments <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
    
    chromosome_list <- report_samp$ref_clust_chr
    correction_factor <- report_samp$correction_factor
    
    
    gg <- ggplot(Granges_segments) +
      ggtitle(paste0("Sample name: ",samples[i]),
              subtitle = paste0("Correction factor= ", correction_factor %>% round(3)," / diploid chrs: ", chromosome_list)) +
      geom_segment(stat="identity", aes(y=CN, colour="old_CN"), size = 1.2, alpha=0.7) +
      geom_segment(stat="identity", aes(y=Granges_segments$CN_correced, colour="new_CN_corrected"), size = 1.2, alpha=0.7) +
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
    
    if(plot_save ==TRUE){
      
      if ( plot_format=="png" ) { png( paste0(plot_path, samples[i],".png"), width = 16, height = 4, units = "in", res = 300 ) 
      } else if(plot_format=="tiff") { tiff( paste0(plot_path, samples[i],".tif"), width = 16, height = 4, units = "in", res = 300 )
      } else if(plot_format=="pdf") { pdf( file = paste0(plot_path, samples[i],".pdf"), width = 16, height = 4 )
      } else if(plot_format == "jpg") { jpeg( paste0(plot_path, samples[i],".jpg"), width = 16, height = 4, units = "in", res = 300 ) 
      } else { message("wrong plot format") 
        break}
      
      print(gg)
      
      
      dev.off()
    }
    
    
    if(plot_viewer==TRUE){
      print(gg)
    }
    
    
  }
}