#' Analyze a data set for gene expression stability with the NormFinder algorithm
#' 
#' @param data A data frame with sample, gene and group (optional) identifier and expression values. 
#' Expression values can be raw Cq or on linear scale.
#' @param group Logical, if TRUE, the function will consider a grouping variable
#' @param sample Character Column name for the sample identifiers
#' @param gene Character Column name of the gene identifiers
#' @param groups Character Column name of the group identifiers
#' @param expression Character Name of column containing expression values 
#' @param log Logical if TRUE the function will consider expression values to be on a log scale 
#' @export
tidy_normfinder<-function(data, 
                     group = TRUE, 
                     sample = "sample", 
                     gene = "gene", 
                     groups = "group", 
                     expression = "expression",
                     log = FALSE){
 
  
  if(log == FALSE) data[, expression] <- log2(data[, expression])
  
  d <- data.frame(sample = data[, sample], 
                  gene = data[, gene], 
                  group = data[, groups],
                  expression = data[, expression])
  

  # Calculate number of groups and number of genes
  if(group == TRUE) {
    if(length(unique(d$group)) < 2) stop("There are less than two groups")
    tbl_group_gene <- table(d$group, d$gene)
    ngroups <- nrow(tbl_group_gene) 
    ngenes <- ncol(tbl_group_gene)
    
    # If groups
    variance_group_gene <- d %>%
      # mean expression per sample
      dplyr:: group_by(sample)%>%
      dplyr::mutate(y_j = mean(expression, na.rm=TRUE))%>%
      # mean expression per gene and group
      dplyr::group_by(group, gene)%>%
      dplyr::mutate(y_ig = mean(expression, na.rm=TRUE))%>%
      # group average
      dplyr:: group_by(group)%>%
      dplyr::mutate(y_g = mean(expression, na.rm=TRUE))%>%
      # mean expression per gene 
      group_by(gene)%>%
      mutate(y_i = mean(expression, na.rm=TRUE))%>%
      # global mean expression
      dplyr::ungroup() %>%
      dplyr::mutate(ybar = mean(expression, na.rm=TRUE))%>%
      dplyr:: group_by(gene, group)%>%
      dplyr::mutate(N = dplyr::n(), 
             r_igj = sum((expression - y_ig - y_j + y_g)^2, na.rm = TRUE)/(N-1)) %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(var.sum.group = sum(unique(r_igj)))%>%
      dplyr::ungroup() %>%
      dplyr::mutate(var = (r_igj - var.sum.group / (ngenes * (ngenes - 1))) / (1 - 2/ngenes)) %>%
      dplyr::group_by(gene, group) %>%
      dplyr::summarise(y_ig = mean(y_ig),
                y_g = mean(y_g),
                y_i = mean(y_i),
                ybar = mean(ybar), 
                N = mean(N), 
                var = mean(var))%>%
      dplyr::ungroup() %>%
      dplyr::mutate(dif = y_ig - y_g - y_i + ybar,
             var_n = var / N,
             sumdd = sum(dif*dif),
             meanvar = mean(var_n),
             n = ((ngroups - 1) * (ngenes - 1)),
             t = max(sumdd / n - meanvar, 0),
             dn = dif * t / (t + var_n),
             varnew = var_n + t * var_n /(t + var_n),
             qm = abs(dn) + sqrt(varnew))%>% 
      dplyr:: group_by(gene)%>%
      dplyr::summarise(rho = mean(qm))%>%
      dplyr::arrange(rho)
    
    return(variance_group_gene)
  }
  if(group == FALSE) {
    
    ngenes <- length(unique(d$gene))
    
    variance_gene <- d %>%
      # mean expression per sample
      dplyr:: group_by(sample)%>%
      dplyr::mutate(y_j = mean(expression, na.rm=TRUE))%>%
      # mean expression per gene 
      dplyr:: group_by(gene)%>%
      dplyr:: mutate(y_i = mean(expression, na.rm=TRUE))%>%
      # global mean expression
      dplyr:: ungroup() %>%
      dplyr:: mutate(ybar = mean(expression, na.rm=TRUE))%>%
      dplyr:: group_by(gene)%>%
      dplyr:: mutate(N = dplyr::n(), 
             r_igj = sum((expression - y_i - y_j + ybar)^2, na.rm = TRUE)/(N-1)) %>%
      dplyr:: ungroup() %>%
      dplyr::mutate(var.sum = sum(unique(r_igj)))%>%
      dplyr::ungroup() %>%
      dplyr::mutate(var = (r_igj - var.sum / (ngenes * (ngenes - 1))) / (1 - 2/ngenes)) %>%
      dplyr:: group_by(gene) %>%
      dplyr::summarise(SDgroup = sqrt(mean(var)))%>%
      dplyr:: arrange(SDgroup)
    
    return(variance_gene)
    
  }
}
  
