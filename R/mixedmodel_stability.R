#' Analyze reference gene stability in a mixed model framework -- Alternative function
#' 
#' The function is an implementation of the method proposed by Dai, H., et al. (2013) with modifications.
#' The function determines gene stability by iterative fitting of potential gene combinations and calculation of a bootstrap confidence interval 
#' around the intraclass correlation. Inference about fixed effects are by default drawn from LRT-tests based where 
#' the full model is compared to a reduced model not containing the effects of interest.
#' 
#' @references Dai, H., et al. (2013). "Mixed modeling and sample size calculations for identifying housekeeping genes." Stat Med 32(18): 3115-3125.
#' @param data a data.frame() in tidy format, each row an observation and each column a variable 
#' @param target character name of column containing gene/target identifier
#' @param response character Specifying the respons column, typically expression values.
#' @param fixed.effects character vector specifying fixed effects included in the model
#' @param random.effect character vector specifying the random effect term, should be specified as used in lme4 e.g. (1|participant)
#' @param interactions Logical, if TRUE, interactions between fixed effects and targets are tested in likelihood ratio tests 
#' @param form formula Specification of formula used in lme4::lmer() for the full model
#' @param reduced.model formula Specification of the model used in calculation of LRT and intraclass correlation
#' @param icc.model formula specification for icc calculation, if NULL then the reduced.model is used
#' @param hypothesis.test Specifies what test to use for the hypothesis test for fixed effects. 
#' The default is likelihood-ratio test ("LRT"), other alternatives are tests of 
#' regression coefficients using the Walds t-test ("wald.t") with Satterthwaite approximation of degrees of 
#' freedom as implemented in the lmerTest package, analysis of variance ("anova") with Satterthwaite approximation 
#' of degrees of freedom as implemented in the lmerTest package or t-test for regression parameters ("z.as.t") specifying a 
#' critical t-value (default 1.96 z- as t-distribution). Type II Wald chisquare tests implemented the car package ("wald.chi").
#' @param LRT.type Character specifying how LRT test should be performed. "global" tests a model containing 
#' all fixed effects against null model (only containing target) and returns a p-value for this test. "add" tests addition 
#' of each fixed effects term to a null-model and returns p-values associated with addition of each term. "interaction" tests
#' two models against reduced models interaction between systematic effects and targets are tested together against a model only 
#' containing fixed effects and fixed effects are tested against a null-model. P-values for interactions and fixed effects are returned.   
#' @param p.threshold numeric Specifying the p-value threshold for hypothesis testing in LRT and Wald t-test of fixed effects.
#' @param critical.t numeric Specifies the t-value threshold for detecting fixed effects in the mixed linear model using 
#' t-test of regression coefficients when hypothesis test is "z.as.t", defaults to ~1.96 or qnorm(0.975).
#' @param icc.interval numeric a single fraction specifying the bootstrap confidence interval, defaults to 0.95 
#' @param icc.type character specifying which method to use for bootstrap confidence interval calculation 
#' ("norm", "basic" or "perc"), default = "norm". See boot::boot.ic for details. 
#' @param n.genes numeric A vector specifying the number of genes in possible combinations of genes to be evaluated. 
#' Can be a vector of e.g. c(2,3), which will perform the algorithm using combinations of two and three genes/targets.
#' @param n.sims numeric Specifies how many bootstraps to perform for the calculation of CI. This process is time consuming when large numbers are used. Defaults to 500 simulations.
#' @param progress Logical default to TRUE gives a progressbar.
#' @return A data frame with with bootstrap confidence intervals for each possible combination of n.genes without significant fixed effects. 
#' @export
mixedmodel_stability<-function(data, 
                               target, 
                               response,
                               fixed.effects,
                               random.effect,
                               form, 
                               reduced.model, 
                               icc.model = NULL, 
                               hypothesis.test = "LRT", 
                               LRT.type = "global",
                               p.threshold = 0.05,
                               critical.t = 1.959964, 
                               icc.interval = 0.95, 
                               icc.type="norm", 
                               n.genes=2, 
                               n.sims=500, 
                               progress=TRUE){
  
  results <- list()
  
  for(n in 1:length(n.genes)){
  
    # possible combinations of size n.genes
    combinations <- combn(unique(data[, which(colnames(data) == target)]), n.genes[n])
    
    # function for calculating ICC based on user defined random effects and bootstrap model
    icc.calc <- function(model){
      icc <- as.numeric(data.frame(VarCorr(model))[1,4]) / sum(as.numeric(data.frame(VarCorr(model))[,4]))
      return(icc)
    }
    
    # Convert to formulas 
    form <- as.formula(form)
    reduced.model <- as.formula(reduced.model)
    if(is.null(icc.model)) icc.model <- as.formula(reduced.model) else icc.model <- as.formula(icc.model)
    
    # A result data frame for icc
    results.icc <- data.frame(gene.combination=rep(NA, ncol(combinations)),
                              icc=rep(NA, ncol(combinations)),
                              icc.l=rep(NA, ncol(combinations)),
                              icc.u=rep(NA, ncol(combinations)))
    
    
    
    # If LRT.type == add
    if(hypothesis.test == "LRT") {
      if(LRT.type == "add") {
        
        response.target <- paste0(response, "~", target)
        systematic <- paste0(fixed.effects, collapse = "+")
        interactions <- paste(target, fixed.effects, sep=":")
        
        terms <- c(fixed.effects, interactions)
        
        results.lrt <- data.frame(matrix(ncol = length(terms), nrow = ncol(combinations)))
        colnames(results.lrt) <- terms 
      }
      
      if(LRT.type == "global") {
        results.lrt <- data.frame(p.val = rep(NA, ncol(combinations)))
      }
      
      if(LRT.type == "interactions") {
        results.lrt <- data.frame(syst.pval = rep(NA, ncol(combinations)),
                                  inter.pval = rep(NA, ncol(combinations)))
      }
      
    }
    
    if(hypothesis.test == "wald.t") {
      
      results.wald.t <- data.frame(coef = rep(NA, ncol(combinations)),
                                   t.val = rep(NA, ncol(combinations)),
                                   p.val = rep(NA, ncol(combinations)))
      
    }
    
    if(hypothesis.test == "z.as.t") {
      
      results.z.as.t <- data.frame(coef = rep(NA, ncol(combinations)),
                                   t.val = rep(NA, ncol(combinations)))
      
    }
    
    

    
    if(progress==TRUE){
      pb <- txtProgressBar(min = 0, max = ncol(combinations), style = 3)
    }
    
    
    for(i in 1:ncol(combinations)){
      
      genes <- combinations[,i] # genes in the specified combination
      results.icc[i,1] <- paste(genes[1:n.genes[n]], collapse=":") # gene combination 
      subset.data <- data[data[, which(colnames(data) == target)] %in% genes,] # subset of genes used in iteration
      
      
      ## Hypothesis testing ##
      if(!(hypothesis.test %in% c("LRT", "wald.t", "z.as.t", "anova", "wald.chi"))) stop("Please specify a hypothesis test") # Do not continue with wrong 

      # 1. LRT 
      
      if(hypothesis.test == "LRT") {
        
        # create formula terms based specifications
        response.target <- paste0(response, "~", target)
        systematic <- paste0(fixed.effects, collapse = "+")
        interactions <- paste(target, fixed.effects, sep=":", collapse = "+")
        
        # combine to formula
        full.f <- as.formula(paste(response.target, systematic, interactions, random.effect, sep = "+"))
        systematic.f <- as.formula(paste(response.target, systematic, random.effect, sep = "+"))
        reduced.f <- as.formula(paste(response.target, random.effect, sep = "+"))
        
        if(LRT.type == "interactions") {
         
          # fit models based on formulas
          full.mod <- lme4::lmer(full.f, data=subset.data, REML=FALSE)
          systematic.mod <- lme4::lmer(systematic.f, data=subset.data, REML=FALSE)
          reduced.mod <- lme4::lmer(reduced.f, data=subset.data, REML=FALSE)
          
          # test for systematic (fixed) effects and their interaction with targets 
          # fixed effects are grouped and tested against a reduced model
          systematic.effects.p <- data.frame(anova(reduced.mod, systematic.mod))[2,8]
          interaction.effects.p <- data.frame(anova(systematic.mod, full.mod))[2,8]
          
          results.lrt[i, 1] <- systematic.effects.p
          results.lrt[i, 2] <- interaction.effects.p
        } 
        
        if(LRT.type == "add") { # test addition of each fixed effect term against a reduced model 
          
          reduced.mod <- lme4::lmer(reduced.f, data = subset.data, REML = FALSE)
          results.lrt[i, ] <-  data.frame(add1(reduced.mod, scope = terms, test = "Chi"))[-1, 4]
          
        } 
        
        if(LRT.type == "global") { # combine all fixed effects in one lrt-test 
          
          reduced.mod <- lme4::lmer(reduced.f, data = subset.data, REML = FALSE)
          full.mod <- lme4::lmer(full.f, data=subset.data, REML=FALSE)
          results.lrt[i, 1] <- data.frame(anova(full.mod, reduced.mod))[2,8]
        
        }
          
          
          if(any(results.lrt[i, ] < p.threshold)) h.test <- TRUE else h.test <- FALSE
        
      }
      
      # 2. Walds t-test
      if(hypothesis.test == "wald.t") {
        
        # response.target <- paste0(response, "~", target)
        # systematic <- paste0(fixed.effects, collapse = "+")
        # interactions <- paste(target, fixed.effects, sep=":", collapse = "+")
        # 
        # full.f <- as.formula(paste(response.target, systematic, interactions, random.effect, sep = "+"))
        
        # Using the full formula to fit a model
        full.model <- lmerTest::lmer(form, data=subset.data, REML=FALSE) 
       # # Regression coefficients are retrieved from summary function
        coef.table <- data.frame(summary(full.model)$coef)
       
        if(ncol(coef.table) == 5) comp.ok <- TRUE else comp.ok <- FALSE
        
        
        if(comp.ok == FALSE) {
          
        warning(paste("Computational difficulties in lmerTest, inference based on critical t-value =", critical.t))
          
          coef.table$coefs <- rownames(coef.table)
          colnames(coef.table) <- c("estimate", "std.error", "t.val", "coefs")
          
          coef.table <- coef.table[!(coef.table$coefs %in% c(paste0(target, genes), "(Intercept)")), ]
          
          results.wald.t[i, 1] <- coef.table[which.max(abs(coef.table$t.val)) , 4] 
          results.wald.t[i, 2] <- coef.table[which.max(abs(coef.table$t.val)) , 3] 
          results.wald.t[i, 3] <- NA
          
          if(any(results.wald.t[i, 2] > critical.t)) h.test <- TRUE else h.test <- FALSE
        }
        
        if(comp.ok == TRUE) { 
        
          coef.table$coefs <- rownames(coef.table)
        
          colnames(coef.table) <- c("estimate", "std.error", "df", "t.val", "p.val", "coefs")
        
        
          coef.table <- coef.table[!(coef.table$coefs %in% c(paste0(target, genes), "(Intercept)")), ]

        
          results.wald.t[i, 1] <- coef.table[which.min(coef.table$p.val) , 6] 
          results.wald.t[i, 2] <- coef.table[which.min(coef.table$p.val) , 4] 
          results.wald.t[i, 3] <- coef.table[which.min(coef.table$p.val) , 5]
          
        # are there significant fixed effects?
        if(results.wald.t[i, 3] < p.threshold) h.test <- TRUE else h.test <- FALSE
        }
      }
      
      # Anova table
      if(hypothesis.test == "anova") {
        
        # Using the full formula to fit a model
        full.model <- lmerTest::lmer(form, data=subset.data, REML=FALSE) 
        # Regression coefficients are retrieved from summary function
        coef.table <- data.frame(anova(full.model))
        
        # Exclude fixed effect of target gene and possibly Intercept (if not )from p-values
        coef.table$coefs <- rownames(coef.table)
        p.vals <- coef.table[!(coef.table$coef %in% c(target)), 6]
        
        # is there a fixed effect
        if(any(p.vals < p.threshold)) h.test <- TRUE else h.test <- FALSE
        
      }
      
      # z as t-distribution for regression coefficients 
      if(hypothesis.test == "z.as.t") {
        
        # create formula terms based specifications
        response.target <- paste0(response, "~", target)
        systematic <- paste0(fixed.effects, collapse = "+")
        interactions <- paste(target, fixed.effects, sep=":", collapse = "+")
        
        # combine to formula
        full.f <- as.formula(paste(response.target, systematic, interactions, random.effect, sep = "+"))
        
        # Using the full formula to fit a model
        full.model <- lme4::lmer(full.f, data=subset.data, REML=TRUE) 
        # Regression coefficients are retrieved from summary function
        coef.table <- data.frame(summary(full.model)$coef)
        
        coef.table$coefs <- rownames(coef.table)
        colnames(coef.table) <- c("estimate", "std.error", "t.val", "coefs")
        
        coef.table <- coef.table[!(coef.table$coefs %in% c(paste0(target, genes), "(Intercept)")), ]
        
        results.z.as.t[i, 1] <- coef.table[which.max(abs(coef.table$t.val)) , 4] 
        results.z.as.t[i, 2] <- coef.table[which.max(abs(coef.table$t.val)) , 3] 
        
        # is there a fixed effect
        if(abs(results.z.as.t[i, 2]) > critical.t) h.test <- TRUE else h.test <- FALSE
        
      }
      # Anova table
      if(hypothesis.test == "wald.chi") {
        
        # Using the full formula to fit a model
        full.model <- lme4::lmer(form, data=subset.data, REML=FALSE) 
        # Analysis of Deviance Table (Type III Wald chisquare tests) from car package
        coef.table <- data.frame(car::Anova(full.model, type = 2))
        
        # Exclude fixed effect of target gene and possibly Intercept (if not )from p-values
        coef.table$coefs <- rownames(coef.table)
        p.vals <- coef.table[!(coef.table$coef %in% c(target, "(Intercept)")), 3]
        
        # is there a fixed effect
        if(any(p.vals < p.threshold)) h.test <- TRUE else h.test <- FALSE
        
      }

      # If an hypothesis test is significant
      # do not calculate bootstrap ICC
      if(h.test == TRUE) {
        results.icc[i,2]<-NA
        results.icc[i,3]<-NA
        results.icc[i,4]<-NA
      } 
      
      else {
        # If no t-values over critical t
        # calculate bootstrap ICC
        
        # fit model for ICC calculations
        icc.mod <- lme4::lmer(icc.model, REML = TRUE, data = subset.data)
        
        b <- lme4::bootMer(icc.mod, icc.calc, use.u = FALSE, nsim = n.sims)
        
        ci <- boot::boot.ci(b, conf = icc.interval, type = icc.type)
        
        if(icc.type=="basic"){
          cis<-c(data.frame(ci[4])[1,4], data.frame(ci[4])[1,5])
        }
        if(icc.type=="norm"){
          cis<-c(data.frame(ci[4])[1,2], data.frame(ci[4])[1,3])
        }
        if(icc.type=="perc"){
          cis<-c(data.frame(ci[4])[1,4], data.frame(ci[4])[1,5])
        }
        
        results.icc[i, 2] <- b$t0
        results.icc[i, 3] <- cis[1]
        results.icc[i, 4] <- cis[2]
        
      }
      
      setTxtProgressBar(pb, i)
    }
    
    if(progress==TRUE) {close(pb)}
    
    if(hypothesis.test == "LRT"){
      r <- cbind(results.icc, results.lrt)
      r$n.genes <- n.genes[n]
      results[[n]] <- r
    } 
    
    if(hypothesis.test == "wald.t") {
      r <- cbind(results.icc, results.wald.t)
      r$n.genes <- n.genes[n]
      results[[n]] <- r
    }
    
    if(hypothesis.test == "z.as.t") {
      r <- cbind(results.icc, results.z.as.t)
      r$n.genes <- n.genes[n]
      results[[n]] <- r
    }
    
    
    if(!(hypothesis.test %in% c("LRT", "wald.t", "z.as.t"))) {
      r <- results.icc
      r$n.genes <- n.genes[n]
      results[[n]] <- r
    }
    
  }
  results <- dplyr::bind_rows(results)
  return(results)
}
