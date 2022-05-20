#' DRrefit_plot
#'
#' @description The function plot the copy number profile before and after DRrefit recalibration 
#'
#' @param corrected_segments DRrefit output dataframe. 
#' @param DRrefit_report  DRrefit output dataframe.
#' @param plot_viewer Logical parameter. When it is TRUE, the function print the output plot in the R viewer.By default is FALSE.
#' @param plot_save  Logical parameter. When it is TRUE, the function save the plot in the chosen path and format. By default is FALSE.
#' @param plot_format File format for the output plots (accepts "png", "jpg", "pdf", "tiff"). By default is "png" 
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
#' data("TCGA_BRCA_CN_segments")
#' 
#' chr_list <- c("10q","11p","12p","19q","1p","21q","2q","3p","4p","4q","6p","6q","7p" )
#' 
#' results <- DRrefit(segments_chort = TCGA_BRCA_CN_segments, chrlist = chr_list)
#'                    
#' my_segments <- results$corrected_segments
#' my_report <- results$report
#' 
#' DRrefit_plot(corrected_segments = my_segments,
#'              DRrefit_report = my_report, 
#'              plot_viewer= FALSE, 
#'              plot_save = FALSE)
#' 
#' 

DRrefit_plot <- function(corrected_segments,
                         DRrefit_report,
                         plot_viewer= F,
                         plot_save = F,
                         plot_format = "png",
                         plot_path
) {
  
  ID <- CN <- NULL
  
  samples <- unique(corrected_segments$ID)
  
  for (i in seq_along(samples)) {
    
    segments <- corrected_segments %>%  filter(ID==samples[i])
    
    report_samp <- DRrefit_report %>%  filter(sample== samples[i])
    
    Granges_segments <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
    
    chromosome_list <- report_samp$ref_clust_chr
    correction_factor <- report_samp$correction_factor
    
    if(abs(correction_factor) > 0.1 ){
      color1= "green4"
      color2="red3"
    } else {
      color1="orange"
      color2="red3"
    } 
    
    gg <- ggplot() +
      ggtitle(paste0("Sample name: ",samples[i]),
              subtitle = paste0("Correction factor= ", correction_factor %>% round(3)," / diploid chrs: ", chromosome_list)) +
      geom_segment(data= Granges_segments,  stat="identity", aes(y= CN),  colour=color1, size = 1.2, alpha=0.7) +
      geom_segment(data = Granges_segments, stat="identity", aes(y= CN_corrected), colour= color2, size = 1.2, alpha=0.7) +
      ylim(0,5) +
      facet_grid( ~seqnames, scales = "free", space = "free", margins = FALSE, switch="x") +
      theme_bw() +
      theme(panel.spacing=unit(.05, "lines"), 
            panel.border = element_rect(color = "black", fill = NA, size = 0.2), 
            axis.text.x = element_text(size=6, angle = 90, vjust = 0.5, hjust=1),
            strip.background = element_rect(color="black", fill="grey90", linetype="solid"),
            strip.text.x = element_text(face = "bold", size=8),
            legend.position="bottom") +
      scale_x_continuous(labels = function(x) format(x/1000000, scientific = FALSE), 
                         breaks = seq(0, 300*10^6,  by = 25*10^6),
                         expand = c(0, 0), 
                         limits = c(NA, NA))+
      geom_hline(yintercept = 2, linetype=3) +
      xlab("Mbp (genomic position)")
    
    
    if(plot_save){
      
      if ( plot_format=="png" ) { png( paste0(plot_path, samples[i],".png"), width = 16, height = 4, units = "in", res = 300 ) 
      } else if(plot_format=="tiff") { tiff( paste0(plot_path, samples[i],".tif"), width = 16, height = 4, units = "in", res = 300 )
      } else if(plot_format=="pdf") { pdf( file = paste0(plot_path, samples[i],".pdf"), width = 16, height = 4 )
      } else if(plot_format == "jpg") { jpeg( paste0(plot_path, samples[i],".jpg"), width = 16, height = 4, units = "in", res = 300 ) 
      } else { message("wrong plot format") 
        break}
      
      print(gg)
      
      
      dev.off()
    }
    
    if(plot_viewer){
      print(gg)
    }
    
    
  }
}