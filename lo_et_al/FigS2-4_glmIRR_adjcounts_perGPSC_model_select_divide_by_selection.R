require(tidyverse)
library(MASS)

#input population size per year from https://github.com/rgladstone/GPSCs/blob/master/lo_et_al/pop_years.csv
pop <- read.csv("pop_years.csv", header = TRUE, sep =",")
#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary.csv", header = TRUE, sep =",")
period <- c("Post-PCV7","Post-PCV13")

#######Calculate IRR per GPSC for pre- vs PCV7 or PCV13 periods############
GPSC_IRR_pre_PCV13 <- matrix(data=NA,nrow=0,ncol=20)
GPSC_IRR_pre_PCV7 <- matrix(data=NA,nrow=0,ncol=20)
for (country in unique(pop$Country)){
  for (post in period){
    #case selection and population size per year for pre-PCV and either post-PCV7 or post-PCV13
    one_pop <- subset(pop, 
                      Country==country & Period==post |
                      Country==country & Period=="Pre-PCV")
    #select rows from country for pre-PCV and either post-PCV7 or post-PCV13
    dat <- subset(paper2,
                    Country==country & Vaccine_Period==post |
                    Country==country & Vaccine_Period=="Pre-PCV", select = c(GPSC, Year))
    GPSCs <- unique(dat$GPSC)
    #Calculate IRR pre-PCV vs post-PCV13 for each country
    for (cluster in GPSCs){
      GPSC_dat <- droplevels(subset(dat, GPSC==cluster, select = c(Year)))
      tab <- as.data.frame(table(GPSC_dat))
      #combine with genome counts per year in NVT-GPSCs/VT-GPSCs
      pop_tab <- merge(one_pop, tab, by.y = "GPSC_dat", by.x = "Year", all.x=TRUE)
      pop_tab$Period <- factor(pop_tab$Period, levels = c("Pre-PCV", post))
      pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
      #adjust genome count by case selection
      pop_tab$estimated.cases <- pop_tab$Freq/pop_tab$selection
      #capture estimated cases per period before rounding
      pre_counts <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Freq)))
      post_counts <- sum(subset(pop_tab, Period==post, select=c(Freq)))
      pre_cases <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(estimated.cases)))
      post_cases <- sum(subset(pop_tab, Period==post, select=c(estimated.cases)))
      pre_years <- dim(subset(pop_tab, Period=="Pre-PCV"))[1]
      post_years <- dim(subset(pop_tab, Period==post))[1]
      pre_population_avg <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(population)))/pre_years
      post_population_avg <- sum(subset(pop_tab, Period==post, select=c(population)))/post_years
      #round estimated counts
      pop_tab$estimated.cases <- round(pop_tab$estimated.cases)
      
      if  (sum(pop_tab$Freq)>=5){
        #Add one estimated case to pre or post if no genomes
        if (sum(subset(pop_tab, Period=="Pre-PCV")['Freq'])==0 | sum(subset(pop_tab, Period==post)['Freq'])==0){
          pop_tab$estimated.cases <- pop_tab$estimated.cases+1
          add <- 1
        } else {
          add <- 0
        }
        
        #calculate IRR using period averages
        IRR_by2 <- matrix(c(post_cases/post_years,pre_cases/pre_years,post_population_avg,pre_population_avg), nrow = 2, byrow = TRUE)
        rownames(IRR_by2) <- c("post_avg_population", "pre_avg_population"); colnames(IRR_by2) <- c("post_est_annual_cases", "pre_est_annual_cases")
        IRR_by2 <- round(IRR_by2)
        res <- epi.2by2(IRR_by2, method = "cross.sectional", conf.level = 0.95, units = 100, homogeneity = "breslow.day",
                        outcome = "as.rows")  
        IRR_calc <- res$res$IRR.strata.wald$est
        confi_lo_calc <- res$res$IRR.strata.wald$lower
        confi_up_calc <- res$res$IRR.strata.wald$upper
        ps_calc <- res$res$chisq.strata$p.value
        
        #runglm
        res = glm(estimated.cases ~ Period + offset(log(population)) , data=pop_tab, family = "poisson")
        #test fit
        GoFit <- 1 - pchisq(summary(res)$deviance, 
                            summary(res)$df.residual)
        #set zero inflation values to NA for models with no zeros
        ZI <- NA
        ZI_period <- NA
        
        #assess zero inflation total and period
        if (sum(pop_tab$estimated.cases ==0)>0){
          res.zip = zeroinfl(estimated.cases ~ Period|Period, data = pop_tab)
          ZI <- summary(res.zip)$coefficients$zero[1,4]
          ZI_period <- summary(res.zip)$coefficients$zero[2,4]
        }

        #If poisson does fit
        if (GoFit >0.05){
          model <- "poisson"
          converged <- res$converged
          pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
          IRR <- unname(pre_post[2])
          ps <-  unname(coef(summary(res))[,4][2])
          pre_postconfint <- try(confint(res) %>% exp)
          confi_lo <- pre_postconfint[2,1]
          confi_up <- pre_postconfint[2,2]
          #average incidence before vaccine (intercept)/100,000 population
          pre_inc <- pre_post[1]*100000
          #average incidence post vaccine/100,000 population derived from pre incidence and IRR
          post_inc <- (pre_post[1]*100000)*pre_post[2]
          
          #If poisson does NOT quite fit (minor violation)
        } else if (GoFit >=0.01 & GoFit <=0.05){
          #Calcultate robust standard errors
          model <- "poisson robust SE"
          converged <- res$converged
          pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
          IRR <- unname(pre_post[2])
          cov.res <- vcovHC(res, type="HC0")
          std.err <- sqrt(diag(cov.res))
          r.est <- cbind(Estimate= coef(res), "Robust SE" = std.err,
                         "Pr(>|z|)" = 2 * pnorm(abs(coef(res)/std.err), lower.tail=FALSE),
                         LL = coef(res) - 1.96 * std.err,
                         UL = coef(res) + 1.96 * std.err)
          ps <- r.est[2,3]
          #convert confidence intervals from log values
          confi_lo <- r.est[2,4] %>% exp
          confi_up <- r.est[2,5] %>% exp
          #average incidence before vaccine (intercept)/100,000 population
          pre_inc <- pre_post[1]*100000
          #average incidence post vaccine/100,000 population derived from pre incidence and IRR
          post_inc <- (pre_post[1]*100000)*pre_post[2]
          
          #If poisson does NOT fit
        } else {
          #fit negative bionomial (this gives theta.ml warning if zero inflation present, which is the next model to be fitted if glm.nb is <0.05)
          model <- "negative bionomial"
          res.nb = glm.nb(estimated.cases ~ Period + offset(log(population)) , data=pop_tab)
          #test fit
          GoFit <- 1 - pchisq(summary(res.nb)$deviance,
                              summary(res.nb)$df.residual)
          converged <- res.nb$converged
          pre_post <- res.nb$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
          IRR <- unname(pre_post[2])
          ps <-  unname(coef(summary(res.nb))[,4][2])
          pre_postconfint <- try(confint(res.nb) %>% exp)
          confi_lo <- pre_postconfint[2,1]
          confi_up <- pre_postconfint[2,2]
          #average incidence before vaccine (intercept)/100,000 population
          pre_inc <- pre_post[1]*100000
          #average incidence post vaccine/100,000 population derived from pre incidence and IRR
          post_inc <- (pre_post[1]*100000)*pre_post[2]
          #assess zero inflation total and period
          if (sum(pop_tab$estimated.cases ==0)>0){
            res.zip = zeroinfl(estimated.cases ~ Period|Period, data = pop_tab)
            ZI <- summary(res.zip)$coefficients$zero[1,4]
            ZI_period <- summary(res.zip)$coefficients$zero[2,4]
          }
          #If negative bionomial does NOT fit
          if (GoFit <0.05 & !is.na(ZI) & ZI <0.05){
            #fit zero inflated negative binomial
            res.zinb = zeroinfl(estimated.cases ~ Period|1 + offset(log(population)) , data=pop_tab, dist="negbin")
            #test fit using Log(theta)
            GoFit <- summary(res.zinb)$coefficients$count[3,4]
            #select model based on theta
            if(summary(res.zinb)$coefficients$count[3,4] <0.05){
              model <- "zero inflated negative binomial"
              converged <- res.zinb$converged
              ZI <- summary(res.zinb)$coefficients$zero[1,4]
              pre_post <- res.zinb$coefficients$count %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
              IRR <- unname(pre_post[2])
              ps <-  unname(coef(summary(res.zinb))$count[,4][2])
              pre_postconfint <- try(confint(res.zinb) %>% exp)
              confi_lo <- pre_postconfint[2,1]
              confi_up <- pre_postconfint[2,2]
              #average incidence before vaccine (intercept)/100,000 population
              pre_inc <- pre_post[1]*100000
              #average incidence post vaccine/100,000 population derived from pre incidence and IRR
              post_inc <- (pre_post[1]*100000)*pre_post[2]
            } else {
              model <- "zero inflated poisson"
              res.zip = zeroinfl(estimated.cases ~ Period|1, data = pop_tab)
              ZI <- summary(res.zip)$coefficients$zero[1,4]
              converged <- res.zip$converged
              pre_post <- res.zip$coefficients$count %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
              IRR <- unname(pre_post[2])
              ps <-  unname(coef(summary(res.zip))$count[,4][2])
              pre_postconfint <- try(confint(res.zip) %>% exp)
              confi_lo <- pre_postconfint[2,1]
              confi_up <- pre_postconfint[2,2]
              #average incidence before vaccine (intercept)/100,000 population
              pre_inc <- pre_post[1]*100000
              #average incidence post vaccine/100,000 population derived from pre incidence and IRR
              post_inc <- (pre_post[1]*100000)*pre_post[2]
           
            }   
          }
        }
        if (post == "Post-PCV7"){
          GPSC_IRR_pre_PCV7 <- rbind(GPSC_IRR_pre_PCV7,c(country, cluster, model, GoFit, ZI, ZI_period, IRR, confi_lo, confi_up, ps, pre_counts, pre_cases, post_counts, post_cases,pre_inc, post_inc,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
          } else {
          GPSC_IRR_pre_PCV13 <- rbind(GPSC_IRR_pre_PCV13,c(country, cluster, model, GoFit, ZI, ZI_period, IRR, confi_lo, confi_up, ps, pre_counts, pre_cases, post_counts, post_cases,pre_inc, post_inc,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
        }
      }
    }
  }
}
#adjust average IRR p-value
GPSC_IRR_pre_PCV7 <- cbind(GPSC_IRR_pre_PCV7,p.adjust(GPSC_IRR_pre_PCV7[,10], method = "BH", n=length(GPSC_IRR_pre_PCV7[,ncol(GPSC_IRR_pre_PCV7)])))
GPSC_IRR_pre_PCV13 <- cbind(GPSC_IRR_pre_PCV13,p.adjust(GPSC_IRR_pre_PCV13[,10], method = "BH", n=length(GPSC_IRR_pre_PCV13[,ncol(GPSC_IRR_pre_PCV13)])))
#adjust model IRR p-value
GPSC_IRR_pre_PCV7 <- cbind(GPSC_IRR_pre_PCV7,p.adjust(GPSC_IRR_pre_PCV7[,20], method = "BH", n=length(GPSC_IRR_pre_PCV7[,ncol(GPSC_IRR_pre_PCV7)])))
GPSC_IRR_pre_PCV13 <- cbind(GPSC_IRR_pre_PCV13,p.adjust(GPSC_IRR_pre_PCV13[,20], method = "BH", n=length(GPSC_IRR_pre_PCV13[,ncol(GPSC_IRR_pre_PCV13)])))
#reorder columns
GPSC_IRR_pre_PCV7 <- cbind(GPSC_IRR_pre_PCV7[,1:10],GPSC_IRR_pre_PCV7[,21],GPSC_IRR_pre_PCV7[,11:20],GPSC_IRR_pre_PCV7[,22])
GPSC_IRR_pre_PCV13 <- cbind(GPSC_IRR_pre_PCV13[,1:10],GPSC_IRR_pre_PCV13[,21],GPSC_IRR_pre_PCV13[,11:20],GPSC_IRR_pre_PCV7[,22])

colnames(GPSC_IRR_pre_PCV7) <- c("Country","GPSC","model","GOF","ZI", "ZI_period","IRR","lower","upper","p","adj.p","pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence","IRR_average","lower","upper","p-value","adj.avg.p")
colnames(GPSC_IRR_pre_PCV13) <- c("Country","GPSC","model","GOF","ZI", "ZI_period","IRR","lower","upper","p","adj.p","pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence","IRR_average","lower","upper","p-value","adj.avg.p")
write.csv(GPSC_IRR_pre_PCV7, file ="FigS2-4_postPCV7_glmIRR_pseudo_adjfreq_model_select_divide_by_selection.csv", row.names = FALSE)
write.csv(GPSC_IRR_pre_PCV13 , file ="FigS2-4_postPCV13_glmIRR_pseudo_adjfreq_model_select_divide_by_selection.csv", row.names = FALSE)
