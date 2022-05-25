#' Popeye
#' @description The function assign the chromosomal arm to each segment.
#'
#' @param segments data.frame formatted with correct column names (see package vignette)
#'
#' @return Return a data frame containg segments with the arm annotation.
#' @export
#' 
#' @importFrom tidyr %>%
#' @importFrom GenomicRanges makeGRangesFromDataFrame seqnames start end
#' @importFrom dplyr filter arrange rename
#' @importFrom plyranges join_overlap_intersect
#' 
#' 
#' @examples
#' data("TCGA_BRCA_CN_segments")
#' data <- TCGA_BRCA_CN_segments[1:9] #as it already presents the arm column
#' data_annotated <- Popeye(segments = data)

Popeye <- function(segments) {
  chr <- ID <- NULL
  
  segments <- segments %>% filter(!chr %in% c("X", "Y"))
  
  chrtab <- BOBaFIT:::chrtab

  chrtabGR <- makeGRangesFromDataFrame(chrtab) 
  
  segmentsGR <- makeGRangesFromDataFrame(segments, keep.extra.columns = TRUE)
  
  segmentsGR_noCentromeres <- join_overlap_intersect(chrtabGR, segmentsGR)

  segments_clean <- segmentsGR_noCentromeres %>% as.data.frame()
  
  arms <- chrtab$chrarm %>% unique()
  
  annot_data <- data.frame()
  i=1
  for( i in seq_along(arms)){
    
    message(i, " - arm: ", arms[i])
    
    chr_sel <- chrtab$chr[i]
    arm_sel <- chrtab$arm[i]
    chrarm_sel <- chrtab$chrarm[i]
    start_sel <- chrtab$start[i]
    end_sel <- chrtab$end[i]
    
    data_f <- segments_clean %>% filter(seqnames==chr_sel, start >= start_sel, end<= end_sel)
    
    
    data_f$arm <- arm_sel
    data_f$chrarm <- chrarm_sel
    
    annot_data <- rbind(annot_data, data_f)
    
  }

  annot_data_sort <- annot_data %>% arrange(ID, seqnames, start)
  
  OUTPUT <- annot_data_sort %>% rename("chr"="seqnames")
  OUTPUT
  
}