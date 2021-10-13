#' computeNormalChromosomes
#'
#' @description This function compute the DRrefits' input "chromosome list". It is a vector that contains the chromosomal arms considered "normal" in the cohort of samples tested (BED file), under a specific tolerance value
#'
#' @param segments data.frame formatted with correct column names
#' @param tolerance_val decimal value of alteration frequency, by default is 0.15
#' @param maxCN threshold of max copy number to consider, by default is 6
#' @param min_threshold minimum threshold to define a normal CN, by default is 1.60
#' @param max_threshold maximum threshold to define a normal CN, by default is 2.40
#'
#' @importFrom dplyr filter group_by summarise arrange
#' @importFrom tidyr %>%
#' @importFrom stats weighted.mean
#' @importFrom ggplot2 aes geom_hline geom_text geom_bar ggplot scale_fill_manual
#' @importFrom stringr str_sort
#'
#'
#' @return vector with chromosome names and plot with the alteration rate of each chromosomal arms
#' @export
#'
#' @examples
#' data(segments)
#' computeNormalChromosomes(segments)



computeNormalChromosomes <- function(segments, tolerance_val = 0.15, maxCN= 6, min_threshold= 1.60, max_threshold=2.40) {

  ID <- arm <- weighted.mean <- CN <- width <- weighted_mean_CN <- alteration_rate <- NULL
  samples <- segments$ID%>% unique()
  segments$CN[segments$CN > 6] <- maxCN

  all_chromosome <- data.frame()

  for (i in seq_along(samples)) {
    message(i," - ",samples[i])
    segments_sample <- segments %>% filter(ID == samples[i])

    CN_CHR <- segments_sample %>%
      group_by(arm) %>%
      summarise(weighted_mean_CN = weighted.mean(CN, w = width),.groups = 'drop_last') %>%
      filter(weighted_mean_CN <= max_threshold & weighted_mean_CN >= min_threshold)

    all_chromosome <- rbind(all_chromosome, CN_CHR)

  }

  result <- (table(all_chromosome$arm)/length(samples) )
  result_filtered <- result[result > 1 - tolerance_val]



  df <- data.frame(arm=names(result), alteration_rate= as.numeric(1- result), observations= table(all_chromosome$arm))
  df$arm <- df$arm %>% factor(levels = str_sort(df$arm, numeric = TRUE))
  print(
    ggplot(df, aes(x=arm, y=alteration_rate)) +
      geom_bar(stat = "identity", aes( fill= alteration_rate > tolerance_val)) +
      geom_hline(yintercept = tolerance_val, colour= "black", linetype= 2)+
      geom_text( aes(label= round(alteration_rate, 2), vjust = 1.4))+
      scale_fill_manual(values = c("#1E90FF","#FF4040")) + 
      theme(legend.position="bottom")
  )

  names(result_filtered)

}
