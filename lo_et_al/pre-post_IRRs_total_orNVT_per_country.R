require(tidyverse)

#input population size per year from https://github.com/rgladstone/GPSCs/blob/master/lo_et_al/pop_years.csv
pop <- read.csv("pop_years.csv", header = TRUE, sep =",")
#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary_v4.csv", header = TRUE, sep =",")

countries <- c("South Africa","USA","Israel")

Total_IRR_PCV13 <- matrix(data=NA,nrow=0,ncol=11)
for (x in countries){
  dat <- subset(paper2, Country==x & Vaccine_Period!="Post-PCV7", select = c(Year))
  tab <- as.data.frame(table(dat))
  one_pop <- subset(pop, 
                    Country==x & Period=="Post-PCV13" |
                      Country==x & Period=="Pre-PCV")
  pop_tab <- merge(one_pop, tab, by.y = "dat", by.x = "Year", all.x=TRUE)
  pop_tab$Period <- factor(pop_tab$Period, levels = c("Pre-PCV", "Post-PCV13"))
  pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
  pop_tab$Actual <- pop_tab$Freq/pop_tab$selection
  #capture estimated cases per period before rounding
  pre_counts <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Freq)))
  post_counts <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(Freq)))
  pre_actual <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Actual)))
  post_actual <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(Actual)))
  #round to estimated counts
  pop_tab$Actual <- round(pop_tab$Actual)
  res = glm(Actual ~ Period + offset(log(population)) , data=pop_tab, family = "poisson")
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
  Total_IRR_PCV13 <- rbind(Total_IRR_PCV13,c(x,IRR,confi_lo,confi_up,ps, pre_counts, pre_actual, post_counts, post_actual,pre_inc, post_inc))
}
colnames(Total_IRR_PCV13) <- c("Country","IRR","lower","upper","p", "pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence")
write.csv(Total_IRR_PCV13, file ="Total_postPCV13_glmIRR.csv", row.names = FALSE)

Total_IRR_PCV13_NVT <- matrix(data=NA,nrow=0,ncol=11)
for (x in countries){
  NVTdat <- subset(paper2, Country==x & Vaccine_Period!="Post-PCV7" & Vaccine_Status=="NVT", select = c(Year))
  tab <- as.data.frame(table(NVTdat))
  one_pop <- subset(pop, 
                    Country==x & Period=="Post-PCV13" |
                      Country==x & Period=="Pre-PCV")
  pop_tab <- merge(one_pop, tab, by.y = "NVTdat", by.x = "Year", all.x=TRUE)
  pop_tab$Period <- factor(pop_tab$Period, levels = c("Pre-PCV", "Post-PCV13"))
  pop_tab$Freq[is.na(pop_tab$Freq)] <- 0
  pop_tab$Actual <- pop_tab$Freq/pop_tab$selection
  #capture estimated cases per period before rounding
  pre_counts <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Freq)))
  post_counts <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(Freq)))
  pre_actual <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Actual)))
  post_actual <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(Actual)))
  #round to estimated counts
  pop_tab$Actual <- round(pop_tab$Actual)
  res = glm(Actual ~ Period + offset(log(population)) , data=pop_tab, family = "poisson")
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
  Total_IRR_PCV13_NVT <- rbind(Total_IRR_PCV13_NVT,c(x,IRR,confi_lo,confi_up,ps, pre_counts, pre_actual, post_counts, post_actual,pre_inc, post_inc))
}
colnames(Total_IRR_PCV13_NVT) <- c("Country","IRR","lower","upper","p", "pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence")
write.csv(Total_IRR_PCV13_NVT, file ="Total_postPCV13_glmNVTIRR.csv", row.names = FALSE)

