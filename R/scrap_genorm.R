
#
#set.seed(1)
#d<-data.frame(sample=seq(1:10),g1= 1/(2^rnorm(10, 23, 1)), 
#              g2=1/(2^rnorm(10, 25,2)), 
#              g3=1/(2^ rnorm(10, 20, 0.1)), 
#              g4=1/(2^rnorm(10, 18, 0.5)))
#
#d <- gather(d, target, expression, g1:g4)
#
#
#library(readxl)
#vand<-read_excel("./data/vandesompele2002.xlsx")
#head(vand)
#
#vand <- gather(vand, gene, expression, ACTB:YWHAZ)
#colnames(vand) <- c("sample", "target", "expression")
#
#tidy_genorm<-function(d, sample, target, expression) {
#  
#  # create relative expression matrix
#  
#  d <- d[, c("sample", "target", "expression")]
#  
#  gene.names <- unique(d[, "target"])
#  
#  d <- tidyr::spread(d, target, expression)
#  
#  dmat <- as.matrix(d[,-1])
#  
#  dmat <- sweep(dmat, 2, apply(dmat, 2, max), FUN="/")
#  
#  pairwise.ratios <- function(x){
#    n <- ncol(x)
#    cn <- colnames(x)
#    
#    cmb <- combn(n, 2)
#    r1 <- apply(cmb, 2, function(j) log(x[, j[1]]/x[, j[2]]))
#    r2 <- apply(cmb, 2, function(j) log(x[, j[2]]/x[, j[1]]))
#    colnames(r1) <- apply(cmb, 2, function(j) paste(cn[j], collapse="."))
#    colnames(r2) <- apply(cmb, 2, function(j) paste(cn[rev(j)], collapse="."))
#    cbind(r1, r2)[, order(c(colnames(r1), colnames(r2)))]
#  }   
#  
#  exclude <- vector()
#  results <- data.frame(gene=rep(NA, length(gene.names)-1), meanM=rep(NA, length(gene.names)-1))
#  
#  
#  
#  
#
#  for(i in 1:length(gene.names)-1){
#    temp<-pairwise.ratios(dmat[, !(colnames(dmat) %in% exclude)])
#    sds <- apply(temp, 2, sd)
#    
#    temp2 <- gather(data.frame(as.list(sds)), gene, sd)
#  
#    temp2$gene <- sapply(strsplit(temp2$gene, "[.]"), `[`, 1)
#
#    m.values <- gather(data.frame(as.list(tapply(temp2$sd, temp2$gene, mean))), gene, m)
#  
#      
#    results[i, 1] <-  m.values[which.max(m.values$m),1]
#    results[i, 2] <- mean(temp2$sd)
#  
#    exclude[i] <-  m.values[which.max(m.values$m),1]
#  }
#  
#  results
#  
#}  
#  
#tidy_genorm(vand[vand$sample %in% c(paste("BM", seq(1:9), sep=""))], "sample", "target", "expression")  
#  
#  vand
#  
#