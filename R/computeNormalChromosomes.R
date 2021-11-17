#' computeNormalChromosomes
#'
#' @description This function compute the DRrefits' input "chromosome list". It is a vector that contains the chromosomal arms considered "normal" in the cohort of samples tested (BED file), under a specific tolerance value
#'
#' @param segments data.frame formatted with correct column names
#' @param tolerance_val decimal value of alteration frequency. By default is 0.15
#' @param maxCN threshold of max copy number to consider. By default is 6
#' @param min_threshold minimum threshold to define a normal CN. By default is 1.60
#' @param max_threshold maximum threshold to define a normal CN. By default is 2.40
#' @param verbose print information about the processes of the function. By default is FALSE
#'
#' @importFrom dplyr filter group_by summarise arrange
#' @importFrom tidyr %>%
#' @importFrom stats weighted.mean
#' @importFrom ggplot2 aes geom_hline geom_text geom_bar ggplot scale_fill_manual labs
#' @importFrom stringr str_sort
#'
#'
#' @return vector with chromosome names and plot with the alteration rate of each chromosomal arms
#' @export
#'
#' @examples
#' data("TCGA_BRCA_CN_segments")
#' chr_list <- computeNormalChromosomes(segments = TCGA_BRCA_CN_segments)



computeNormalChromosomes <- function(segments,
                                     tolerance_val = 0.15,
                                     maxCN = 6,
                                     min_threshold = 1.60,
                                     max_threshold = 2.40,
                                     verbose = FALSE) {
  ID <-
    arm <-
    chrarm <-
    weighted.mean <-
    CN <- width <- weighted_mean_CN <- alteration_rate <- NULL
  
  samples <- segments$ID %>% unique()
  segments$CN[segments$CN > 6] <- maxCN
  
  all_chromosome <- data.frame()
  
  if (verbose == TRUE) {  
    for (i in seq_along(samples)) {
      message(i, " - ", samples[i])
      segments_sample <- segments %>% filter(ID == samples[i])
      
      CN_CHR <- segments_sample %>%
        group_by(chrarm) %>%
        summarise(weighted_mean_CN = weighted.mean(CN, w = width),
                  .groups = 'drop_last') %>%
        filter(weighted_mean_CN <= max_threshold &
                 weighted_mean_CN >= min_threshold)
      
      all_chromosome <- rbind(all_chromosome, CN_CHR)
      
    }
  } else {
    for (i in seq_along(samples)) {
      segments_sample <- segments %>% filter(ID == samples[i])
      
      CN_CHR <- segments_sample %>%
        group_by(chrarm) %>%
        summarise(weighted_mean_CN = weighted.mean(CN, w = width),
                  .groups = 'drop_last') %>%
        filter(weighted_mean_CN <= max_threshold &
                 weighted_mean_CN >= min_threshold)
      
      all_chromosome <- rbind(all_chromosome, CN_CHR)
    }
  }
  
  result <- (table(all_chromosome$chrarm) / length(samples))
  result_filtered <- result[result > 1 - tolerance_val]
  
  
  
  df <-
    data.frame(
      arm = names(result),
      alteration_rate = as.numeric(1 - result),
      observations = table(all_chromosome$chrarm)
    )
  df$arm <-
    df$arm %>% factor(levels = str_sort(df$arm, numeric = TRUE))
  print(
    ggplot(df, aes(x = arm, y = alteration_rate)) +
      geom_bar(stat = "identity", aes(fill = alteration_rate > tolerance_val)) +
      geom_hline(
        yintercept = tolerance_val,
        colour = "black",
        linetype = 2
      ) +
      geom_text(aes(
        label = round(alteration_rate, 2), vjust = 1.4
      )) +
      scale_fill_manual(values = c("#1E90FF", "#FF4040")) +
      theme(legend.position = "bottom") +
      labs(y = "Alteration Rate (%)", x = "Chromosomal arm")
  )
  
  names(result_filtered)
  
}
