require(tidyverse)
library(MASS)
library(pscl)
library(sandwich)

#input population size per year from https://github.com/rgladstone/GPSCs/blob/master/lo_et_al/pop_years.csv
pop <- read.csv("pop_years.csv", header = TRUE, sep =",")
#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary.csv", header = TRUE, sep =",")

countries <- c("South Africa","USA","Israel")

Total_IRR_PCV13 <- matrix(data=NA,nrow=0,ncol=18)
for (x in countries){
  dat <- subset(paper2, Country==x & Vaccine_Period!="Post-PCV7", select = c(Year))
  tab <- as.data.frame(table(dat))
  one_pop <- subset(pop, 
                    Country==x & Period=="Post-PCV13" |
                      Country==x & Period=="Pre-PCV")
  pop_tab <- merge(one_pop, tab, by.y = "dat", by.x = "Year", all.x=TRUE)
  pop_tab$Period <- factor(pop_tab$Period, levels = c("Pre-PCV", "Post-PCV13"))
  pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
  pop_tab$estimated.cases <- pop_tab$Freq/pop_tab$selection
  #capture estimated cases per period before rounding
  pre_counts <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Freq)))
  post_counts <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(Freq)))
  pre_cases <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(estimated.cases)))
  post_cases <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(estimated.cases)))
  pre_years <- dim(subset(pop_tab, Period=="Pre-PCV"))[1]
  post_years <- dim(subset(pop_tab, Period=="Post-PCV13"))[1]
  pre_population_avg <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(population)))/pre_years
  post_population_avg <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(population)))/post_years
  #round to estimated counts
  pop_tab$estimated.cases <- round(pop_tab$estimated.cases)
 
  #2by2 IRR using period averages
  IRR_by2 <- matrix(c(post_cases/post_years,pre_cases/pre_years,post_population_avg,pre_population_avg), nrow = 2, byrow = TRUE)
  rownames(IRR_by2) <- c("post_avg_population", "pre_avg_population"); colnames(IRR_by2) <- c("post_est_annual_cases", "pre_est_annual_cases")
  IRR_by2 <- round(IRR_by2)
  res <- epi.2by2(IRR_by2, method = "cross.sectional", conf.level = 0.95, units = 100, homogeneity = "breslow.day",
                  outcome = "as.rows")  
  IRR_calc <- res$res$IRR.strata.wald$est
  confi_lo_calc <- res$res$IRR.strata.wald$lower
  confi_up_calc <- res$res$IRR.strata.wald$upper
  ps_calc <- res$res$chisq.strata$p.value
  
  #run poisson
  res = glm(estimated.cases ~ Period + offset(log(population)) , data=pop_tab, family = "poisson")
  #test fit
  GoFit <- 1 - pchisq(summary(res)$deviance, 
             summary(res)$df.residual)
  if (GoFit >0.05){
    model <- "poisson"
    converge <- res$converged
    pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
    IRR <- unname(pre_post[2])
    ps <-  unname(coef(summary(res))[,4][2])
    pre_postconfint <- try(confint(res) %>% exp)
    confi_lo <- pre_postconfint[2,1]
    confi_up <- pre_postconfint[2,2]
  } else if (GoFit >=0.01 & GoFit <=0.05){
    #Calcultate robust standard errors
    converge <- res$converged
    pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
    IRR <- unname(pre_post[2])
    cov.res <- vcovHC(res, type="HC0")
    std.err <- sqrt(diag(cov.res))
    r.est <- cbind(Estimate= coef(res), "Robust SE" = std.err,
                   "Pr(>|z|)" = 2 * pnorm(abs(coef(res)/std.err), lower.tail=FALSE),
                   LL = coef(res) - 1.96 * std.err,
                   UL = coef(res) + 1.96 * std.err)
    ps <- r.est[2,3]
    #calculate convert confidence intervals from log values
    confi_lo <- r.est[2,4] %>% exp
    confi_up <- r.est[2,5] %>% exp
    model <- "poisson robust SE"
  } else {
    #fit negative bionomial
    res.nb = glm.nb(estimated.cases ~ Period + offset(log(population)), data=pop_tab)
    #test fit
    GoFit <- 1 - pchisq(summary(res.nb)$deviance,
               summary(res.nb)$df.residual)
    model <- "negative bionomial"
    converge <- res.nb$converged
    pre_post <- res.nb$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
    IRR <- unname(pre_post[2])
    ps <-  unname(coef(summary(res.nb))[,4][2])
    pre_postconfint <- try(confint(res.nb) %>% exp)
    confi_lo <- pre_postconfint[2,1]
    confi_up <- pre_postconfint[2,2]
    }

  #average incidence before vaccine (intercept)/100,000 population
  pre_inc <- pre_post[1]*100000
  #average incidence post vaccine/100,000 population derived from pre incidence and IRR
  post_inc <- (pre_post[1]*100000)*pre_post[2]
  Total_IRR_PCV13 <- rbind(Total_IRR_PCV13,c(x,model,converge,GoFit,IRR,confi_lo,confi_up,ps, pre_counts, pre_cases, post_counts, post_cases,pre_inc, post_inc,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
}
colnames(Total_IRR_PCV13) <- c("Country","model","Converged","GOF", "IRR","lower","upper","p-value","pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence","IRR_average","lower","upper","p-value")
write.csv(Total_IRR_PCV13, file ="Total_postPCV13_glmIRR_model_select_divide_by_select.csv", row.names = FALSE)

Total_IRR_PCV13_NVT <- matrix(data=NA,nrow=0,ncol=18)
for (x in countries){
  NVTdat <- subset(paper2, Country==x & Vaccine_Period!="Post-PCV7" & Vaccine_Status=="NVT", select = c(Year))
  tab <- as.data.frame(table(NVTdat))
  one_pop <- subset(pop, 
                    Country==x & Period=="Post-PCV13" |
                      Country==x & Period=="Pre-PCV")
  pop_tab <- merge(one_pop, tab, by.y = "NVTdat", by.x = "Year", all.x=TRUE)
  pop_tab$Period <- factor(pop_tab$Period, levels = c("Pre-PCV", "Post-PCV13"))
  pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
  pop_tab$estimated.cases <- pop_tab$Freq/pop_tab$selection
  #capture estimated cases per period before rounding
  pre_counts <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Freq)))
  post_counts <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(Freq)))
  pre_cases <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(estimated.cases)))
  post_cases <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(estimated.cases)))
  pre_years <- dim(subset(pop_tab, Period=="Pre-PCV"))[1]
  post_years <- dim(subset(pop_tab, Period=="Post-PCV13"))[1]
  pre_population_avg <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(population)))/pre_years
  post_population_avg <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(population)))/post_years
  #round to estimated counts
  pop_tab$estimated.cases <- round(pop_tab$estimated.cases)
  
  #2by2 IRR using period averages
  IRR_by2 <- matrix(c(post_cases/post_years,pre_cases/pre_years,post_population_avg,pre_population_avg), nrow = 2, byrow = TRUE)
  rownames(IRR_by2) <- c("post_avg_population", "pre_avg_population"); colnames(IRR_by2) <- c("post_est_annual_cases", "pre_est_annual_cases")
  IRR_by2 <- round(IRR_by2)
  res <- epi.2by2(IRR_by2, method = "cross.sectional", conf.level = 0.95, units = 100, homogeneity = "breslow.day",
                  outcome = "as.rows")  
  IRR_calc <- res$res$IRR.strata.wald$est
  confi_lo_calc <- res$res$IRR.strata.wald$lower
  confi_up_calc <- res$res$IRR.strata.wald$upper
  ps_calc <- res$res$chisq.strata$p.value
  
  #run poisson
  res = glm(estimated.cases ~ Period + offset(log(population)) , data=pop_tab, family = "poisson")
  #test fit
  GoFit <- 1 - pchisq(summary(res)$deviance, 
                      summary(res)$df.residual)
  if (GoFit >0.05){
    model <- "poisson"
    converge <- res$converged
    pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
    IRR <- unname(pre_post[2])
    ps <-  unname(coef(summary(res))[,4][2])
    pre_postconfint <- try(confint(res) %>% exp)
    confi_lo <- pre_postconfint[2,1]
    confi_up <- pre_postconfint[2,2]
  } else if (GoFit >=0.01 & GoFit <=0.05){
    #Calcultate robust standard errors
    model <- "poisson robust SE"
    converge <- res$converged
    pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
    IRR <- unname(pre_post[2])
    cov.res <- vcovHC(res, type="HC0")
    std.err <- sqrt(diag(cov.res))
    r.est <- cbind(Estimate= coef(res), "Robust SE" = std.err,
                   "Pr(>|z|)" = 2 * pnorm(abs(coef(res)/std.err), lower.tail=FALSE),
                   LL = coef(res) - 1.96 * std.err,
                   UL = coef(res) + 1.96 * std.err)
    ps <- r.est[2,3]
    #calculate convert confidence intervals from log values
    confi_lo <- r.est[2,4] %>% exp
    confi_up <- r.est[2,5] %>% exp
  } else {
    #fit negative bionomial
    res.nb = glm.nb(estimated.cases ~ Period + offset(log(population)) , data=pop_tab)
    #test fit
    GoFit <- 1 - pchisq(summary(res.nb)$deviance,
                        summary(res.nb)$df.residual)
    model <- "negative bionomial"
    converge <- res.nb$converged
    pre_post <- res.nb$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
    IRR <- unname(pre_post[2])
    ps <-  unname(coef(summary(res.nb))[,4][2])
    pre_postconfint <- try(confint(res.nb) %>% exp)
    confi_lo <- pre_postconfint[2,1]
    confi_up <- pre_postconfint[2,2]
  }
  #average incidence before vaccine (intercept)/100,000 population
  pre_inc <- pre_post[1]*100000
  #average incidence post vaccine/100,000 population derived from pre incidence and IRR
  post_inc <- (pre_post[1]*100000)*pre_post[2]
  Total_IRR_PCV13_NVT <- rbind(Total_IRR_PCV13_NVT,c(x,model,converge,GoFit,IRR,confi_lo,confi_up,ps,pre_counts, pre_cases, post_counts, post_cases,pre_inc, post_inc,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
}
colnames(Total_IRR_PCV13_NVT) <- c("Country","model","Converged", "GOF","IRR","lower","upper","p-value","pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence","IRR_average","lower","upper","p-value")
write.csv(Total_IRR_PCV13_NVT, file ="TotalNVT_postPCV13_glmIRR_model_select_divide_by_select.csv", row.names = FALSE)

