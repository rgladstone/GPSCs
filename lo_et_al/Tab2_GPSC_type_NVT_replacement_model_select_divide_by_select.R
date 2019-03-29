require(tidyverse)
library(MASS)
library(pscl)
library(sandwich)
library("epiR")

#input population size per year from https://github.com/rgladstone/GPSCs/blob/master/lo_et_al/pop_years.csv
pop <- read.csv("pop_years.csv", header = TRUE, sep =",")
#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary.csv", header = TRUE, sep =",")

#Calculate IRR for NVT in NVT-GPSCs or VT-GPSCs between pre-PCV and post-PCV13 periods for each country
GPSCtype_NVTIRR_pre_PCV13 <- matrix(data=NA,nrow=0,ncol=21)

for (country in unique(pop$Country)){
  #Determine GPSC type per country based on VT proportion (<50% or >=50%) in first period observed
  GPSC_props <- matrix(data=NA,nrow=0,ncol=3)
  colnames(GPSC_props) <- c("GPSC","Proportion VT", "type")
  coun_p2 <- subset(paper2, Country==country, select=c(GPSC,Vaccine_Period,Vaccine_Status,Year))
  total_NVTpre <- dim(subset(coun_p2, Vaccine_Period=="Pre-PCV" & Vaccine_Status=="NVT", select=c("Vaccine_Status")))[1]
  total_NVTpost <- dim(subset(coun_p2, Vaccine_Period=="Post-PCV13" & Vaccine_Status=="NVT", select=c("Vaccine_Status")))[1]
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
  
  #test NVT in NVT-GPSCs then in VT-GPSCs
  for (type in GPSC_type){
    #what type is being tested?
    if ((all.equal(GPSC_type$`NVT-GPSC`, type) =="TRUE")[1]){
      what_type <- "NVT-GPSCs"
    } else {
      what_type <- "VT-GPSCs"
    }
    #How many GPSCs?
    No_GPSCs <- length(type)
        #select year of isolation for NVTs collected pre-PCV7 or post-PCV13 in NVT-GPSCs/VT-GPSCs from country of interest
    NVTdat <- droplevels(subset(paper2, GPSC %in% type & Country==country & Vaccine_Period!="Post-PCV7" & Vaccine_Status=="NVT", select = c(Year)))
    tab <- as.data.frame(table(NVTdat))
    #case selection and population size per year
    one_pop <- subset(pop, 
                      Country==country & Period=="Post-PCV13" |
                        Country==country & Period=="Pre-PCV")
    #combine with genome counts per year in NVT-GPSCs/VT-GPSCs
    pop_tab <- merge(one_pop, tab, by.y = "NVTdat", by.x = "Year", all.x=TRUE)
    pop_tab$Period <- factor(pop_tab$Period, levels = c("Pre-PCV", "Post-PCV13"))
    pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
    #adjust genome count by case selection
    pop_tab$estimated.cases <- pop_tab$Freq/pop_tab$selection
    #capture estimated cases per period before rounding
    pre_counts <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Freq)))
    post_counts <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(Freq)))
    #Proportion of NVTs found in NVT-GPSCs or VT-GPSCs
    prop_of_NVTpre <- (pre_counts/(total_NVTpre))*100
    prop_of_NVTpost <- (post_counts/(total_NVTpost))*100
    pre_cases <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(estimated.cases)))
    post_cases <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(estimated.cases)))  
    pre_years<- dim(subset(pop_tab, Period=="Pre-PCV"))[1]
    post_years<- dim(subset(pop_tab, Period=="Post-PCV13"))[1]
    pre_population_avg <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(population)))/pre_years
    post_population_avg <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(population)))/post_years
    #round estimated counts
    pop_tab$estimated.cases <- round(pop_tab$estimated.cases)
    
    #run poisson
    res = glm(estimated.cases ~ Period + offset(log(population)) , data=pop_tab, family = "poisson")
    #test fit
    GoFit <- 1 - pchisq(summary(res)$deviance, 
                        summary(res)$df.residual)
  
    ZI <- NA
    ZI_period <- NA
    
    #assess zero inflation total and period
    if (sum(subset(pop_tab, Period=="Pre-PCV")$estimated.cases ==0)>0){
      res.zip.period = zeroinfl(estimated.cases ~ Period|Period, data = pop_tab)
      ZI <- summary(res.zip.period)$coefficients$zero[1,4]
      ZI_period <- summary(res.zip.period)$coefficients$zero[2,4]
    }
    if (GoFit >0.05){
      model <- "poisson"
      converged <- res$converged
      pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
      IRR <- unname(pre_post[2])
      ps <-  unname(coef(summary(res))[,4][2])
      pre_postconfint <- try(confint(res) %>% exp)
      confi_lo <- pre_postconfint[2,1]
      confi_up <- pre_postconfint[2,2]
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
      #calculate convert confidence intervals from log values
      confi_lo <- r.est[2,4] %>% exp
      confi_up <- r.est[2,5] %>% exp
      model <- "poisson robust SE"
    } else {
      #fit negative bionomial
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
      #assess zero inflation total and period
      if (sum(subset(pop_tab, Period=="Pre-PCV")$estimated.cases ==0)>0){
        res.zip.period = zeroinfl(estimated.cases ~ Period|Period, data = pop_tab)
        ZI <- summary(res.zip.period)$coefficients$zero[1,4]
        ZI_period <- summary(res.zip.period)$coefficients$zero[2,4]
      }
    }
    if (GoFit<0.04 & model != "poisson robust SE"){
      model <- "none"
      dat <- matrix(c(post_cases/post_years,pre_cases/pre_years,post_population_avg,pre_population_avg), nrow = 2, byrow = TRUE)
      rownames(dat) <- c("post_avg_population", "pre_avg_population"); colnames(dat) <- c("post_est_annual_cases", "pre_est_annual_cases")
      dat <- round(dat)
      res <- epi.2by2(dat = as.table(dat), method = "cross.sectional", conf.level = 0.95, units = 100, homogeneity = "breslow.day",
               outcome = "as.rows")  
      IRR <- res$res$IRR.strata.wald$est
      confi_lo <- res$res$IRR.strata.wald$lower
      confi_up <- res$res$IRR.strata.wald$upper
      ps <- res$res$chisq.strata$p.value
      GoFit <- "NA"
    }
    #average incidence before vaccine (intercept)/100,000 population
    pre_inc <- pre_post[1]*100000
    #average incidence post vaccine/100,000 population derived from pre incidence and IRR
    post_inc <- (pre_post[1]*100000)*pre_post[2]
    #Replacement
    extra_cases <- post_inc-pre_inc
    GPSCtype_NVTIRR_pre_PCV13 <- rbind(GPSCtype_NVTIRR_pre_PCV13,c(country, what_type, No_GPSCs, prop_of_NVTpre, prop_of_NVTpost, model, converged, GoFit, ZI, ZI_period, IRR, confi_lo, confi_up, ps, pre_counts, pre_cases, post_counts, post_cases,pre_inc, post_inc,extra_cases ))
  }
}
colnames(GPSCtype_NVTIRR_pre_PCV13) <- c("Country","GPSC-type", "No. GPSCs","Prop of NVTs pre-PCV","Prop of NVTs post-PCV","model","Converged","GOF","Zero inflation (ZI)", "ZI by period","IRR","lower","upper","p-value", "pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence","extra cases")
write.csv(GPSCtype_NVTIRR_pre_PCV13, file ="Tab2_postPCV13_glmIRR_NVT_GPSCtype_model_select_divide_by_select.csv", row.names = FALSE)
