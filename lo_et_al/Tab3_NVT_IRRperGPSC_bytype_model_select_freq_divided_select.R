require(tidyverse)
library(MASS)
library(pscl)
library(sandwich)

#input population size per year from https://github.com/rgladstone/GPSCs/blob/master/lo_et_al/pop_years.csv
pop <- read.csv("pop_years.csv", header = TRUE, sep =",")
#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary.csv", header = TRUE, sep =",")

#######Calculate NVT IRR per GPSC############
GPSC_NVTIRR_pre_PCV13 <- matrix(data=NA,nrow=0,ncol=25)
GPSC_NVTIRR_pre_PCV7 <- matrix(data=NA,nrow=0,ncol=25)
GPSC_NVTIRR_PCV7_PCV13 <- matrix(data=NA,nrow=0,ncol=25)

for (country in unique(pop$Country)){
  #Determine GPSC type per country based on VT proportion (<50% or >=50%) in first period observed
  GPSC_props <- matrix(data=NA,nrow=0,ncol=3)
  colnames(GPSC_props) <- c("GPSC","Proportion VT", "type")
  coun_p2 <- subset(paper2, Country==country, select=c(GPSC,Vaccine_Period,Vaccine_Status,In_Silico_Serotype,Country,Year))
  GPSCs <- unique(coun_p2$GPSC)
  for (cluster in GPSCs){
    GPSC_data <- subset(coun_p2, GPSC==cluster)
    #table of vaccine period by vaccine status, remove rows that sum to sero (no isolates from that period)
    GPSC_tab <- as.data.frame(unclass(table(GPSC_data$Vaccine_Period,GPSC_data$Vaccine_Status)))[ rowSums(as.data.frame(unclass(table(GPSC_data$Vaccine_Period,GPSC_data$Vaccine_Status))))!=0, ]
    #Calcuate the proportion VT in the earliest period found in last row of table
    GPSC_prop <- tail(GPSC_tab$PCV/rowSums(GPSC_tab), n=1)
    GPSC_props <- rbind(GPSC_props, c(cluster,GPSC_prop,ifelse(GPSC_prop<0.5,"NVT","VT")))
  }
  GPSC_props <- as.data.frame(GPSC_props)
  `NVT-GPSCs` <- as.character(subset(GPSC_props, type=="NVT", select=c(GPSC))$GPSC)
  `VT-GPSCs` <- as.character(subset(GPSC_props, type=="VT", select=c(GPSC))$GPSC)
  GPSC_type <- list(`NVT-GPSCs`,`VT-GPSCs`)
  names(GPSC_type) <- c("NVT-GPSC","VT-GPSC")
  
  #Calculate NVT IRR pre-PCV vs post-PCV13 for each country per GPSC
  period <- c("Post-PCV7","Post-PCV13")
  for (post in period){
    #case selection and population size per year for pre-PCV and either post-PCV7 or post-PCV13
    one_pop <- subset(pop, 
                      Country==country & Period==post |
                        Country==country & Period=="Pre-PCV")
    #select NVT rows from country for pre-PCV and either post-PCV7 or post-PCV13
    dat <- subset(paper2, 
                  Country==country & Vaccine_Period==post & Vaccine_Status=="NVT"|
                    Country==country & Vaccine_Period=="Pre-PCV" & Vaccine_Status=="NVT",
                  select = c(GPSC, Year))
    GPSCs <- unique(dat$GPSC)
    #test IRR for NVT in each GPSC
    for (cluster in GPSCs){
      #what type of GPSC is being tested?
      if (cluster %in% GPSC_type$`NVT-GPSC`){
        what_type <- "NVT-GPSC"
      } else {
        what_type <- "VT-GPSC"
      }
      GPSC_dat <- droplevels(subset(dat, GPSC==cluster, select = c(Year)))
      tab <- as.data.frame(table(GPSC_dat))
      #combine with genome counts per year
      pop_tab <- merge(one_pop, tab, by.y = "GPSC_dat", by.x = "Year", all.x=TRUE)
      pop_tab$Period <- factor(pop_tab$Period, levels = c("Pre-PCV", post))
      pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
      pre_freq <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Freq)))
      post_freq <- sum(subset(pop_tab, Period==post, select=c(Freq)))
      #estimate cases
      pop_tab$estimated.cases <- round(pop_tab$Freq/pop_tab$selection)
      pre_cases <-  sum(subset(pop_tab, Period=="Pre-PCV")['estimated.cases'])
      post_cases <-  sum(subset(pop_tab, Period==post)['estimated.cases'])
      pre_years <- dim(subset(pop_tab, Period=="Pre-PCV"))[1]
      post_years <- dim(subset(pop_tab, Period==post))[1]
      pre_zeros <- sum(subset(pop_tab, Period=="Pre-PCV")$estimated.cases == 0) 
      post_zeros <- sum(subset(pop_tab, Period==post)$estimated.cases == 0)
      pre_population_avg <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(population)))/pre_years
      post_population_avg <- sum(subset(pop_tab, Period==post, select=c(population)))/post_years
      #Only test GPSCs with at least 5  NVT genomes
      if  (sum(pop_tab$Freq)>=5){
        #Add 1 to all case estimate if no genomes observed in one or other period
        if (sum(subset(pop_tab, Period=="Pre-PCV")['Freq'])==0 | sum(subset(pop_tab, Period==post)['Freq'])==0){
          pop_tab$estimated.cases <- pop_tab$estimated.cases+1
          add <- 1
        } else {
          add <- 0
        }
        
        #capture cases after 1 added to all if one period=0
        pre_cases_SN <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(estimated.cases)))
        post_cases_SN <- sum(subset(pop_tab, Period==post, select=c(estimated.cases)))
        
        #calculate IRR using period averages
        IRR_by2 <- matrix(c(post_cases_SN/post_years,pre_cases_SN/pre_years,post_population_avg,pre_population_avg), nrow = 2, byrow = TRUE)
        #calculate IRR using period averages
          rownames(IRR_by2) <- c("post_avg_population", "pre_avg_population"); colnames(IRR_by2) <- c("post_est_annual_cases", "pre_est_annual_cases")
          IRR_by2 <- round(IRR_by2)
          res <- epi.2by2(IRR_by2, method = "cross.sectional", conf.level = 0.95, units = 100, homogeneity = "breslow.day",
                          outcome = "as.rows")  
          IRR_calc <- res$res$IRR.strata.wald$est
          confi_lo_calc <- res$res$IRR.strata.wald$lower
          confi_up_calc <- res$res$IRR.strata.wald$upper
          ps_calc <- res$res$chisq.strata$p.value
        #run model selection routine
        res = glm(estimated.cases ~ Period + offset(log(population)) , data=pop_tab, family = "poisson")
        #test fit
        GoFit <- 1 - pchisq(summary(res)$deviance, 
                            summary(res)$df.residual)
        #set zero inflation to NA as models with no zeros will not be tested for zero inflation
        ZI <- NA
        ZI_period <- NA
        
        #assess zero inflation total and period if any zeros observed
        if (sum(subset(pop_tab)$estimated.cases ==0)>0){
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
          ps <-  r.est[2,3]
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
          if (sum(pop_tab$estimated.cases ==0)>0) {
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
          GPSC_NVTIRR_pre_PCV7 <- rbind(GPSC_NVTIRR_pre_PCV7,c(country, cluster, what_type, pre_years, post_years, pre_zeros, post_zeros, pre_cases, post_cases, pre_inc, post_inc, add, model, converged, GoFit, ZI, ZI_period,IRR, confi_lo, confi_up, ps,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
        } else {
          GPSC_NVTIRR_pre_PCV13 <- rbind(GPSC_NVTIRR_pre_PCV13,c(country, cluster, what_type, pre_years, post_years, pre_zeros, post_zeros,pre_cases, post_cases, pre_inc, post_inc, add,model, converged, GoFit, ZI, ZI_period,IRR, confi_lo, confi_up, ps,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
        }
      }
    }
  }
}





#Calculate NVT IRR post-PCV7 vs post-PCV13 for each country per GPSC
post <- "Post-PCV13"

for (country in unique(pop$Country)){
  #Determine GPSC type per country based on VT proportion (<50% or >=50%) in first period observed
  GPSC_props <- matrix(data=NA,nrow=0,ncol=3)
  colnames(GPSC_props) <- c("GPSC","Proportion VT", "type")
  coun_p2 <- subset(paper2, Country==country, select=c(GPSC,Vaccine_Period,Vaccine_Status,Year))
  GPSCs <- unique(coun_p2$GPSC)
  for (cluster in GPSCs){
    GPSC_data <- subset(coun_p2, GPSC==cluster)
    #table of vaccine period by vaccine status, remove rows that sum to sero (non isolates from that period)
    GPSC_tab <- as.data.frame(unclass(table(GPSC_data$Vaccine_Period,GPSC_data$Vaccine_Status)))[ rowSums(as.data.frame(unclass(table(GPSC_data$Vaccine_Period,GPSC_data$Vaccine_Status))))!=0, ]
    #Calcuate the proportion VT in the earliest period found in last row of table
    GPSC_prop <- tail(GPSC_tab$PCV/rowSums(GPSC_tab), n=1)
    GPSC_props <- rbind(GPSC_props, c(cluster,GPSC_prop,ifelse(GPSC_prop<0.5,"NVT","VT")))
  }
  GPSC_props <- as.data.frame(GPSC_props)
  `NVT-GPSCs` <- as.character(subset(GPSC_props, type=="NVT", select=c(GPSC))$GPSC)
  `VT-GPSCs` <- as.character(subset(GPSC_props, type=="VT", select=c(GPSC))$GPSC)
  GPSC_type <- list(`NVT-GPSCs`,`VT-GPSCs`)
  names(GPSC_type) <- c("NVT-GPSC","VT-GPSC")
  
  post <- "Post-PCV13"
  one_pop <- subset(pop, 
                    Country==country & Period==post |
                      Country==country & Period=="Post-PCV7")
  dat <- subset(paper2, 
                Country==country & Vaccine_Period==post & Vaccine_Status=="NVT"|
                  Country==country & Vaccine_Period=="Post-PCV7"& Vaccine_Status=="NVT",
                select = c(GPSC, Year))
  GPSCs <- unique(dat$GPSC)
  for (cluster in GPSCs){
    #what type of GPSC is being tested?
    if (cluster %in% GPSC_type$`NVT-GPSC`){
      what_type <- "NVT-GPSC"
    } else {
      what_type <- "VT-GPSC"
    }
    
    GPSC_dat <- droplevels(subset(dat, GPSC==cluster, select = c(Year)))
    tab <- as.data.frame(table(GPSC_dat))
    pop_tab <- merge(one_pop, tab, by.y = "GPSC_dat", by.x = "Year", all.x=TRUE)
    pop_tab$Period <- factor(pop_tab$Period, levels = c("Post-PCV7", post))
    pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
    #Estimate cases
    pop_tab$estimated.cases <- round(pop_tab$Freq/pop_tab$selection)
    pre_cases <-  sum(subset(pop_tab, Period=="Post-PCV7")['estimated.cases'])
    post_cases <-  sum(subset(pop_tab, Period==post)['estimated.cases'])
    pre_years <- dim(subset(pop_tab, Period=="Post-PCV7"))[1]
    post_years <- dim(subset(pop_tab, Period==post))[1]
    post7_zeros <- sum(subset(pop_tab, Period=="Post-PCV7")$estimated.cases == 0) 
    post13_zeros <- sum(subset(pop_tab, Period=="Post-PCV13")$estimated.cases == 0) 
    pre_population_avg <- sum(subset(pop_tab, Period=="Post-PCV7", select=c(population)))/pre_years
    post_population_avg <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(population)))/post_years
    #Only test GPSCs with at least 5 NVT genomes
    if  (sum(pop_tab$Freq)>=5){
      #Add 1 to all case estimate if no genomes observed in one or other period
      if (sum(subset(pop_tab, Period=="Post-PCV7")['Freq'])==0 | sum(subset(pop_tab, Period=="Post-PCV13")['Freq'])==0){
        pop_tab$estimated.cases <- pop_tab$estimated.cases+1
        add <- 1
      } else {
        add <- 0
      }
      
      #capture cases after 1 added to all if one period=0
      pre_cases_SN <- sum(subset(pop_tab, Period=="Post-PCV7", select=c(estimated.cases)))
      post_cases_SN <- sum(subset(pop_tab, Period==post, select=c(estimated.cases)))
      
      #calculate IRR using period averages
      IRR_by2 <- matrix(c(post_cases_SN/post_years,pre_cases_SN/pre_years,post_population_avg,pre_population_avg), nrow = 2, byrow = TRUE)
      rownames(IRR_by2) <- c("post_avg_population", "pre_avg_population"); colnames(IRR_by2) <- c("post_est_annual_cases", "pre_est_annual_cases")
      IRR_by2 <- round(IRR_by2)
      res <- epi.2by2(IRR_by2, method = "cross.sectional", conf.level = 0.95, units = 100, homogeneity = "breslow.day",
                      outcome = "as.rows")  
      IRR_calc <- res$res$IRR.strata.wald$est
      confi_lo_calc <- res$res$IRR.strata.wald$lower
      confi_up_calc <- res$res$IRR.strata.wald$upper
      ps_calc <- res$res$chisq.strata$p.value
      
      #run model selection routine
      res = glm(estimated.cases ~ Period + offset(log(population)) , data=pop_tab, family = "poisson")
      #test fit
      GoFit <- 1 - pchisq(summary(res)$deviance, 
                          summary(res)$df.residual)

      #set zero inflation to NA as models with no zeros will not be tested for zero inflation
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
        ps <-  unname(coef(summary(res))[,4][2])
        cov.res <- vcovHC(res, type="HC0")
        std.err <- sqrt(diag(cov.res))
        r.est <- cbind(Estimate= coef(res), "Robust SE" = std.err,
                       "Pr(>|z|)" = 2 * pnorm(abs(coef(res)/std.err), lower.tail=FALSE),
                       LL = coef(res) - 1.96 * std.err,
                       UL = coef(res) + 1.96 * std.err)
        ps <-  r.est[2,3]
        #convert confidence intervals from log values
        confi_lo <- r.est[2,4] %>% exp
        confi_up <- r.est[2,5] %>% exp
        #average incidence before vaccine (intercept)/100,000 population
        pre_inc <- pre_post[1]*100000
        #average incidence post vaccine/100,000 population derived from pre incidence and IRR
        post_inc <- (pre_post[1]*100000)*pre_post[2]
        #assess zero inflation total and period
        if (sum(subset(pop_tab, Period=="Post-PCV7")$estimated.cases ==0)>0 | sum(subset(pop_tab, Period==post)$estimated.cases ==0)>0){
          res.zip = zeroinfl(estimated.cases ~ Period|Period, data = pop_tab)
          ZI <- summary(res.zip)$coefficients$zero[1,4]
          ZI_period <- summary(res.zip)$coefficients$zero[2,4]
        }
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
        if (sum(subset(pop_tab, Period=="Pre-PCV")$estimated.cases ==0)>0){
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
    GPSC_NVTIRR_PCV7_PCV13 <- rbind(GPSC_NVTIRR_PCV7_PCV13,c(country, cluster, what_type, post7_years, post13_years, post7_zeros, post13_zeros, pre_cases, post_cases, pre_inc, post_inc,add, model,converged,GoFit, ZI, ZI_period,IRR, confi_lo, confi_up, ps,IRR_calc,confi_lo_calc,confi_up_calc,ps_calc))
    }
  }
}  


#adjust average IRR p value
GPSC_NVTIRR_pre_PCV7 <- cbind(GPSC_NVTIRR_pre_PCV7,p.adjust(GPSC_NVTIRR_pre_PCV7[,25], method = "BH", n=length(GPSC_NVTIRR_pre_PCV7[,1])))
GPSC_NVTIRR_pre_PCV13 <- cbind(GPSC_NVTIRR_pre_PCV13,p.adjust(GPSC_NVTIRR_pre_PCV13[,25], method = "BH", n=length(GPSC_NVTIRR_pre_PCV13[,1])))
GPSC_NVTIRR_PCV7_PCV13 <- cbind(GPSC_NVTIRR_PCV7_PCV13,p.adjust(GPSC_NVTIRR_PCV7_PCV13[,25], method = "BH", n=length(GPSC_NVTIRR_PCV7_PCV13[,1])))
#adjust average IRR p value
GPSC_NVTIRR_pre_PCV7 <- cbind(GPSC_NVTIRR_pre_PCV7,p.adjust(GPSC_NVTIRR_pre_PCV7[,21], method = "BH", n=length(GPSC_NVTIRR_pre_PCV7[,1])))
GPSC_NVTIRR_pre_PCV13 <- cbind(GPSC_NVTIRR_pre_PCV13,p.adjust(GPSC_NVTIRR_pre_PCV13[,21], method = "BH", n=length(GPSC_NVTIRR_pre_PCV13[,1])))
GPSC_NVTIRR_PCV7_PCV13 <- cbind(GPSC_NVTIRR_PCV7_PCV13,p.adjust(GPSC_NVTIRR_PCV7_PCV13[,21], method = "BH", n=length(GPSC_NVTIRR_PCV7_PCV13[,1])))
colnames(GPSC_NVTIRR_pre_PCV7) <- c("Country","GPSC","GPSC_type","pre_years", "post7_years", "pre_zeros", "post7_zeros","Pre-PCV cases", "Post cases", "pre_incidence", "postPCV7_incidence","add","model", "Converged","GOF","ZI","ZI_period","IRR","lower","upper", "p", "IRR_average","lower","upper","p-value","adj.avg.p","adj.p")
colnames(GPSC_NVTIRR_pre_PCV13) <- c("Country","GPSC","GPSC_type", "pre_years", "post13_years", "pre_zeros", "post13_zeros","Pre-PCV cases","Post cases","pre_incidence","postPCV13_incidence", "add","model", "Converged","GOF","ZI","ZI_period", "IRR","lower","upper", "p","IRR_average","lower","upper","p-value", "adj.avg.p","adj.p")
colnames(GPSC_NVTIRR_PCV7_PCV13) <- c("Country","GPSC","GPSC_type","post7_years", "post13_years", "post7_zeros", "post13_zeros","Pre-PCV cases", "Post cases","postPCV7_inc", "postPCV13_incidence","add","model", "Converged","GOF", "ZI","ZI_period", "IRR","lower","upper", "p","IRR_average","lower","upper","p-value", "adj.avg.p","adj.p")
#reorder cols
GPSC_NVTIRR_pre_PCV7 <- GPSC_NVTIRR_pre_PCV7[,c(1:21,27,22:26)]
GPSC_NVTIRR_pre_PCV13 <- GPSC_NVTIRR_pre_PCV13[,c(1:21,27,22:26)]
GPSC_NVTIRR_PCV7_PCV13 <- GPSC_NVTIRR_PCV7_PCV13[,c(1:21,27,22:26)]

write.csv(GPSC_NVTIRR_pre_PCV7, file ="postPCV7_glmIRR_NVT_pseudo_adjfreq_model_select_divide_by_selection.csv", row.names = FALSE)
write.csv(GPSC_NVTIRR_pre_PCV13 , file ="postPCV13_glmIRR_NVT_pseudo_adjfreq_model_select_divide_by_selection.csv", row.names = FALSE)
write.csv(GPSC_NVTIRR_PCV7_PCV13 , file ="postPCV7_postPCV13_glmIRR_NVT_pseudo_adjfreq_model_select_divide_by_selection.csv", row.names = FALSE)

pre_v_PCV7_PCV7_v_PCV13 <- merge(GPSC_NVTIRR_pre_PCV7[,c(1:3,10:11,15,18:22)],GPSC_NVTIRR_PCV7_PCV13[,c(1,2,10:11,15,18:22)], by=c("Country","GPSC"), all=TRUE)
pre_v_PCV7__PCV7_v_PCV13_pre_v_PCV13 <- merge(pre_v_PCV7_PCV7_v_PCV13,GPSC_NVTIRR_pre_PCV13[,c(1,2,10:11,15,18:22)], by=c("Country","GPSC"), all=TRUE)
colnames(pre_v_PCV7__PCV7_v_PCV13_pre_v_PCV13) <- c("Country","GPSC","GPSC_type","pre_incidence","postPCV7_inc","GOF","IRR-pre_v_PCV7","lower","upper","p","adj.p","postPCV7_inc","postPCV13_incidence","GOF", "IRR-PCV7_v_PCV13","lower","upper","p","adj.p","pre_incidence", "postPCV13_incidence","GOF","IRR-pre_v_PCV13","lower","upper","p","adj.p")
write.csv(pre_v_PCV7__PCV7_v_PCV13_pre_v_PCV13, file ="Tab3_multiperiod_summary_glmIRR_NVT_pseudo_adjfreq_divide_by_selection.csv", row.names = FALSE)

