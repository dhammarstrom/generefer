 #' 
 #' Calculate reference gene stability
 #' 
 #' @description Calculates M values for n >= 2 genes. When there are more than two genes in the data set, 
 #' the function calculates M for the full data set, removes the least stable gene (largest M) and 
 #' calculates M again until only two genes remains. The remaining two genes cannot be ranked. 
 #' Based on gene ranking normalisation factors (NF) are calculated based on two to n genes and pairwise
 #' variations between NFn and NFn+1 are calculated which can be used to determine the optimal
 #' number of reference genes to include in the experiment. 
 #' 
 #' @references Vandesompele, J., et al. (2002). "Accurate normalization of real-time quantitative 
 #' RT-PCR data by geometric averaging of multiple internal control genes." Genome Biol 3(7)  
 #' @param dat Name of the data frame, supplied to the function in i tidy way. 
 #' A column for sample id ("sample"), a column for gene id ("gene") and a column 
 #' containing expression values ("expression").
 #' @param sample Character column name of column containing sample id's
 #' @param gene Character column name of column containing gene id's
 #' @param expression Character column name of column containing sample expression values
 #' @param log Logical, is data on the log scale?
 #' @return A list containing: M, a data frame with stability measures, average stability measures per gene after step-wise
 #' exclusion of genes; Normalisation factors for n >= 2 genes per sample, and; Pairwise variations of normalisation factors. 
 #' @import "dplyr"
 #' @import "tidyr"
 #' @export
 tidy_genorm <- function(dat, sample = "sample", gene = "gene", expression = "expression", log = FALSE) {
   exp.matrix <- dat[ , c(sample, gene, expression)]  # reduce data set to only containing relevant data
   
   exp.matrix <- exp.matrix %>% 
     spread(gene, expression) # spread data
   
   sample.id <- exp.matrix[,1] # store sample ids
   
   exp.matrix <- data.matrix(exp.matrix[,-1]) # remove sample id
   
   n.full <- ncol(exp.matrix) # Calculate n genes
   
   if(n.full < 2) stop(cat(paste("The function only identifies", n.full, "gene. Is your data in the correct format?")))
   
   remove <- vector("character", length = n.full)
   M <- data.frame(gene = rep(NA, length = n.full), 
                   M = rep(NA, length = n.full), 
                   M.avg = rep(NA, length = n.full),
                   n.genes = rep(NA, length = n.full))
   
   if(n.full > 2){   
     for(k in 1:(n.full-2)) { 
       
       
       reduced.ex <- exp.matrix[, !(colnames(exp.matrix) %in% remove)]
       n <- ncol(reduced.ex)
       
       m <- vector(mode = "numeric", length = n)
       names(m) <- colnames(reduced.ex)
       
       
       
       for(i in 1:n) { 
         
         if(log == TRUE){
           r <- reduced.ex[ , i] - reduced.ex[ , -i] # calculates pair-wise ratios stores in r matrix
         } else {
           r <- log2(reduced.ex[ , i] / reduced.ex[ , -i]) # calculates pair-wise ratios stores in r matrix
         }
         
         
         m[i] <-  sum(apply(r, 2, sd)) / (n-1)
         
       }
       
       M[k, 1] <- names(which.max(m)) # Name of gene with least stable 
       M[k, 2] <-  m[which.max(m)] # M value for excluded gene
       M[k, 3] <- mean(m, na.rm = T) # Average M for genes in iteration
       M[k, 4] <- sum(remove == "") # n genes remaining
       
       
       remove[k]<- names(which.max(m))
       
     }}
   
   # Retrieve last pair (most stable)
   last.pair <- exp.matrix[, !(colnames(exp.matrix) %in% remove)]
   
   # Set non-ranked genes as last in the data frame
   M[c(n.full-1, n.full) , 1] <-  colnames(exp.matrix[, !(colnames(exp.matrix) %in% remove)])
   
   # N genes remaining when last pair is calculated
   M[c(n.full-1, n.full) , 4] <- 2
   
   # Calculate SD of pairwise variation in last gene pair
   M[c(n.full-1, n.full), 2] <- sd(log2(last.pair[,1] / last.pair[,2]), na.rm = TRUE) 
   
   M[c(n.full-1, n.full), 3] <- sd(log2(last.pair[,1] / last.pair[,2]), na.rm = TRUE) 
   
   
   genorm_results <- list()
   
   
   nf <- matrix(ncol = n.full - 1, nrow = length(sample.id))
   
   for(n in 0:(n.full-2)) {
     # Calculates geometric average of n reference genes
     nf[,n+1] <- exp(rowMeans(log(exp.matrix[, M[!(M$n.genes > (n.full-n)), 1]]), na.rm = TRUE, dims = 1))
   }
   
   colnames(nf) <- paste0("NF", c(n.full:2))
   
   # reverse order of matrix
   nf <- nf[, ncol(nf):1]
   
   # calculate pairwise variations
   V <- data.frame(NFpair = paste0("V", c(2: (ncol(nf))),"/", c(2: (ncol(nf)) + 1)),
                   pairwise.variation = rep(NA, length = ncol(nf) - 1))
   
   for(i in 1:(ncol(nf)-1)){
     V[i, 2] <- sd(log2(nf[, i] / nf[, i+1]), na.rm=T) # sd of log2 ratios
   }
   
   # combine sample ids with normalisation factors 
   nf <- cbind(data.frame(sample.id = sample.id), data.frame(nf))
   
   ## Results list
   genorm_results[[1]] <- M
   genorm_results[[2]] <- nf
   genorm_results[[3]] <- V
   
   names(genorm_results) <- c("M" , "normalisation.factors", "pairwise.variations")
   
   return(genorm_results)
   
 } # End function
 
 
 
 