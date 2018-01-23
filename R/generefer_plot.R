#' Plot of a generefer data frame
#' 
#' @param dat A data.frame() from lmer_stability()
#' @param icc.threshold Exclude observations based on intra-class correlation lower than threshold, defaults to 0.5
#' @return A list containing a plot and a tidy sorted data frame
#' @import "ggplot2"
#' @import "dplyr"
#' @import "tidyr"
#' @export

generefer_plot <- function(dat, icc.threshold = 0.5) {
 
return.dat  <- dat %>%
    filter(icc > icc.threshold)%>%
    mutate(gene.combination=factor(gene.combination, levels=gene.combination[order(icc.l)]) )

return.plot <- ggplot(return.dat, aes(gene.combination, icc))+
  geom_point()+
  geom_errorbar(aes(ymin=icc.l, ymax=icc.u))+
  coord_flip()
 


return(list(data = return.dat, plot = return.plot))

}


