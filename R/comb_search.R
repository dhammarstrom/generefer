#' Search for gene combinations of size *n_genes* without significant effects of factor variables
#'
#' @param data A data frame with subject,timepoint,condition,target and expression. Expression values should be log transformed linearized Cq values.
#' @param subject numeric column id for subject identifier
#' @param timepoint numeric column id for timepoint identifier
#' @param condition numeric column identifier for condition identifier
#' @param target numeric column identifier for target/gene identifier
#' @param expression numeric column identifier for expression data
#' @param n_genes numeric, specifies how many genes should be combined
#' @import "dplyr"
#' @import "tidyr"
#' @import "lme4"
#' @export
gene_comb<-function(data, subject=1, timepoint=2, condition=3, target=4, expression=5, n_genes=2){
  
  combinations<-combn(unique(data[,target]), n_genes)
  
  potential.gene.comb<-list()
  discarderd.gene.comb<-list()
  models<-list()
  combs<-data.frame(combs=rep(NA, ncol(combinations)))
  
  names(data)[subject]<-"subj"
  names(data)[timepoint]<-"time"
  names(data)[condition]<-"cond"
  names(data)[target]<-"gene"
  names(data)[expression]<-"expr"
  
  
  ## Initialize a Progress Bar
  pb <- txtProgressBar(min = 0, max = ncol(combinations), style = 3) #text based bar
  
  
  for(i in 1:ncol(combinations)){
    
    temp<-data.frame(data[data$gene %in% combinations[,i],])
    
    
    m1<-lme4::lmer(expr~1+(1+gene|subj:cond), data=temp, REML=FALSE)
    m2<-lme4::lmer(expr~time+(1+gene|subj:cond), data=temp, REML=FALSE)
    m3<-lme4::lmer(expr~time+cond+(1+gene|subj:cond), data=temp, REML=FALSE)
    m4<-lme4::lmer(expr~time+cond+(time*cond)+(1+gene|subj:cond), data=temp, REML=FALSE)
    
    
    setTxtProgressBar(pb, i)
    
    if(any(data.frame(anova(m1,m2, m3, m4))[c(2:4),8]<0.05)){discarderd.gene.comb[[i]]<-combinations[,i]}
    
    else {
      combs[i,1]<-i
    }
  }
  
  close(pb)
  
  return(combinations[,!is.na(combs[1])])
}
