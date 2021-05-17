#' computeNormalChromosomes
#'
#' @param segments data.frame formatted with correct column names
#' @param perc decimal value of alteration frequency, by default is 0.15
#' @param maxCN threshold of max copy number to consider, by default is 6
#' @param min_threshold minimum threshold to define a normal CN, by default is 1.60
#' @param max_threshold maximum threshold to define a normal CN, by default is 2.40
#'
#' @importFrom dplyr filter group_by summarise arrange
#' @importFrom tidyr %>%
#' @importFrom stats weighted.mean
#' @importFrom ggplot2 aes geom_hline geom_text geom_bar ggplot
#' @importFrom stringr str_sort
#'
#'
#' @return vector with chormosome names and plot with the alteration rate of each chromosomal arms
#' @export
#'
#' @examples
#' \donttest{
#' data(segments)
#' computeNormalChromosomes(segments)
#' }


computeNormalChromosomes <- function(segments, perc = 0.15, maxCN= 6, min_threshold= 1.60, max_threshold=2.40) {

  ID <- arm <- weighted.mean <- CN <- width <- weighted_mean_CN <- normal_frequency <- NULL
  samples <- segments$ID%>% unique()
  segments$CN[segments$CN > 6] <- maxCN

  segments$width <- segments$end - segments$start

  all_chromosome <- data.frame()

  for (i in 1:length(samples)) {
    print(paste0( i," - ",samples[i]))
    segments_sample <- segments %>% filter(ID == samples[i])

    CN_CHR <- segments_sample %>%
      group_by(arm) %>%
      summarise(weighted_mean_CN = weighted.mean(CN, w = width),.groups = 'drop_last') %>%
      filter(weighted_mean_CN <= max_threshold & weighted_mean_CN >= min_threshold)

    all_chromosome <- rbind(all_chromosome, CN_CHR)

  }

  result <- (table(all_chromosome$arm)/length(samples) )
  result_filtered <- result[result > 1 - perc]



  df <- data.frame(arm=names(result), normal_frequency= as.numeric(1- result), observations= table(all_chromosome$arm))
  df$arm <- df$arm %>% factor(levels = str_sort(df$arm, numeric = TRUE))
  print(
    ggplot(df, aes(x=arm, y=normal_frequency)) +
      geom_bar(stat = "identity", aes( fill= normal_frequency > perc)) +
      geom_hline(yintercept = perc, colour= "black", linetype= 2)+
      geom_text( aes(label= round(normal_frequency, 2), vjust = 1.4))
  )

  names(result_filtered)

}
