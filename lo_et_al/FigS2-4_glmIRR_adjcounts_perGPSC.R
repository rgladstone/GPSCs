#Significant changes measure with IRR via poisson regression within GPSC for South Africa, USA and Israel 

require(tidyverse)

#input population size per year from https://github.com/rgladstone/GPSCs/blob/master/lo_et_al/pop_years.csv
pop <- read.csv("pop_years.csv", header = TRUE, sep =",")
#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary_v4.csv", header = TRUE, sep =",")
period <- c("Post-PCV7","Post-PCV13")

#######Calculate IRR per GPSC for pre- vs PCV7 or PCV13 periods############
GPSC_IRR_pre_PCV13 <- matrix(data=NA,nrow=0,ncol=12)
GPSC_IRR_pre_PCV7 <- matrix(data=NA,nrow=0,ncol=12)
for (country in unique(pop$Country)){
  for (post in period){
    #case selection and population size per year for pre-PCV and either post-PCV7 or post-PCV13
    one_pop <- subset(pop, 
                      Country==country & Period==post |
                      Country==country & Period=="Pre-PCV")
    #select rows from country for pre-PCV and either post-PCV7 or post-PCV13
    dat <- subset(paper2,
                    Country==country & Vaccine_Period==post |
                    Country==country & Vaccine_Period=="Pre-PCV", select = c(GPSC_new, Year))
    GPSCs <- unique(dat$GPSC_new)
    #Calculate IRR pre-PCV vs post-PCV13 for each country
    for (GPSC in GPSCs){
      GPSC_dat <- droplevels(subset(dat, GPSC_new==GPSC, select = c(Year)))
      tab <- as.data.frame(table(GPSC_dat))
      #combine with genome counts per year in NVT-GPSCs/VT-GPSCs
      pop_tab <- merge(one_pop, tab, by.y = "GPSC_dat", by.x = "Year", all.x=TRUE)
      pop_tab$Period <- factor(pop_tab$Period, levels = c("Pre-PCV", post))
      pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
      #adjust genome count by case selection
      pop_tab$Actual <- pop_tab$Freq/pop_tab$selection
      #capture estimated cases per period before rounding
      pre_counts <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Freq)))
      post_counts <- sum(subset(pop_tab, Period==post, select=c(Freq)))
      pre_actual <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Actual)))
      post_actual <- sum(subset(pop_tab, Period==post, select=c(Actual)))
      #round estimated counts
      pop_tab$Actual <- round(pop_tab$Actual)
      if  (sum(pop_tab$Actual)>10){
        if (sum(subset(pop_tab, Period==post)['Actual'])==0){
          pop_tab$Actual <- pop_tab$Actual+1
        }
        if (sum(subset(pop_tab, Period=="Pre-PCV")['Actual'])==0){
          pop_tab$Actual <- pop_tab$Actual+1
        }
        res = glm(Actual ~ Period + offset(log(population)) , data=pop_tab, family = "poisson")
        pre_post <- res$coefficients %>% exp #gives average incidence before vaccine (intercept) & IRR (Post)
        IRR <- unname(pre_post[2])
        ps <-  unname(coef(summary(res))[,4][2])
        pre_postconfint <- try(confint(res) %>% exp)
        if (class(pre_postconfint)=="try-error"){
          confi_lo <- "NA"
          confi_up <- "NA"
          } else {
          confi_lo <- pre_postconfint[2,1]
          confi_up <- pre_postconfint[2,2] 
          }
        #average incidence before vaccine (intercept)/100,000 population
        pre_inc <- pre_post[1]*100000
        #average incidence post vaccine/100,000 population derived from pre incidence and IRR
        post_inc <- (pre_post[1]*100000)*pre_post[2]
        if (post == "Post-PCV7"){
          GPSC_IRR_pre_PCV7 <- rbind(GPSC_IRR_pre_PCV7,c(country, GPSC, IRR, confi_lo, confi_up, ps, pre_counts, pre_actual, post_counts, post_actual,pre_inc, post_inc))
          } else {
          GPSC_IRR_pre_PCV13 <- rbind(GPSC_IRR_pre_PCV13,c(country, GPSC, IRR, confi_lo, confi_up, ps, pre_counts, pre_actual, post_counts, post_actual,pre_inc, post_inc))
        }
      }
    }
  }
}
GPSC_IRR_pre_PCV7 <- cbind(GPSC_IRR_pre_PCV7,p.adjust(GPSC_IRR_pre_PCV7[,6], method = "BH", n=length(GPSC_IRR_pre_PCV7[,6])))
GPSC_IRR_pre_PCV13 <- cbind(GPSC_IRR_pre_PCV13,p.adjust(GPSC_IRR_pre_PCV13[,6], method = "BH", n=length(GPSC_IRR_pre_PCV13[,6])))
#reorder columns
GPSC_IRR_pre_PCV7 <- cbind(GPSC_IRR_pre_PCV7[,1:6],GPSC_IRR_pre_PCV7[,13],GPSC_IRR_pre_PCV7[,7:12])
GPSC_IRR_pre_PCV13 <- cbind(GPSC_IRR_pre_PCV13[,1:6],GPSC_IRR_pre_PCV13[,13],GPSC_IRR_pre_PCV13[,7:12])

colnames(GPSC_IRR_pre_PCV7) <- c("Country","GPSC", "IRR","lower","upper","p","adj.p","pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence")
colnames(GPSC_IRR_pre_PCV13) <- c("Country","GPSC","IRR","lower","upper","p","adj.p","pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence")
write.csv(GPSC_IRR_pre_PCV7, file ="FigS7postPCV7_glmIRR_pseudo_adjfreq.csv", row.names = FALSE)
write.csv(GPSC_IRR_pre_PCV13 , file ="postPCV13_glmIRR_pseudo_adjfreq.csv", row.names = FALSE)
