#' DRrefit_plot
#'
#' @param segments_refitted DRrefit output dataframe 
#' @param DRrefit_report  DRrefit output dataframe 
#' @param plot_format file format for the output plots. By default is "png" (accepts "png", "jpg", "pdf", "tif")
#' @param plot_path path to save output plots
#'
#' @return
#' @export
#' 
#' @importFrom grDevices png dev.off
#' @importFrom ggplot2 ggplot ggtitle ylim facet_grid theme_bw theme scale_x_continuous geom_hline  xlab scale_color_manual unit element_rect
#' @importFrom ggbio geom_segment
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tidyr %>%
#' @importFrom dplyr filter 
#'
#' @examples
#' #' chr_list <- c("10q","11p","12p","19q","1p","21q","2q","3p","4p","4q","6p","6q","7p" )
#' results <- DRrefit(segments,chrlist = chr_list)
#' DRrefit_report <- results$report
#' segments_refitted <- results$segments_corrected
#' DRrefit_plot(segments_refitted,DRrefit_report,plot_path)
#' 
#' 

DRrefit_plot <- function(segments_refitted,
                         DRrefit_report,
                         plot_format = "png",
                         plot_path) {
  
  samples <- unique(segments_refitted$ID)
  i=1
  for (i in seq_along(samples)) {
    
    segments <- segments_refitted %>%  filter(ID==samples[i])
    
    report_samp <- DRrefit_report %>%  filter(sample== samples[i])
    
    Granges_segments <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
    
    chromosome_list <- report_samp$ref_clust_chr
    correction_factor <- report_samp$correction_factor
    
    if ( plot_format=="png" ) { png( paste0(plot_path, samples[i],"_",clust_method,".png"), width = 16, height = 4, units = "in", res = 300 ) 
    } else if(plot_format=="tiff") { tiff( paste0(plot_path, samples[i],"_",clust_method,".tif"), width = 16, height = 4, units = "in", res = 300 )
    } else if(plot_format=="pdf") { pdf( file = paste0(plot_path, samples[i],"_",clust_method,".pdf"), width = 16, height = 4 )
    } else if(plot_format == "jpg") { jpeg( paste0(plot_path, samples[i],"_",clust_method,".jpg"), width = 16, height = 4, units = "in", res = 300 ) 
    } else { message("wrong plot format") 
      break}
    
    print(
      ggplot(Granges_segments) +
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
    )
    
    options(ggplot2. ="viridis")
    
    dev.off()
  }
}