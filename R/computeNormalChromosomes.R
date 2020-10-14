#' computeNormalChromosomes
#'
#' @param segments data.frame formatted with correct column names
#' @param perc decimal value of alteration frequency, by default is 0.15
#'
#' @importFrom dplyr filter group_by summarise arrange
#' @importFrom tidyr %>%
#' @importFrom stats weighted.mean
#'
#'
#' @return vector with chormosome names
#' @export
#'
#' @examples
#' \dontrun{
#' computeNormalChromosomes(segments)
#' }


computeNormalChromosomes <- function(segments, perc = 0.15) {

  ID <- arm <- weighted.mean <- CN <- width <- weighted_mean_CN <- NULL
  samples <- segments$ID%>% unique()

  segments$width <- segments$end - segments$start

  all_chromosome <- data.frame()

  i=1
  for (i in 1:length(samples)) {

    segments_sample <- segments %>% filter(ID == samples[i])

    CN_CHR <- segments_sample %>%
      group_by(arm) %>%
      summarise(weighted_mean_CN = weighted.mean(CN, w = width),.groups = 'drop_last') %>%
      filter(weighted_mean_CN <= 2.10 & weighted_mean_CN >=1.90)

    all_chromosome <- rbind(all_chromosome, CN_CHR)

  }

  result <- (table(all_chromosome$arm)/length(samples) )
  result_filtered <- result[result > 1 - perc]
  names(result_filtered)
}
