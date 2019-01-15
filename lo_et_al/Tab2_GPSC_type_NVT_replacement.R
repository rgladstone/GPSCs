require(tidyverse)

#input population size per year from https://github.com/rgladstone/GPSCs/blob/master/lo_et_al/pop_years.csv
pop <- read.csv("pop_years.csv", header = TRUE, sep =",")
#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary.csv", header = TRUE, sep =",")

#Calculate IRR for NVT in NVT-GPSCs or VT-GPSCs between pre-PCV and post-PCV13 periods for each country
GPSCtype_NVTIRR_pre_PCV13 <- matrix(data=NA,nrow=0,ncol=16)

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
    pop_tab$Actual <- pop_tab$Freq/pop_tab$selection
    #capture estimated cases per period before rounding
    pre_counts <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Freq)))
    post_counts <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(Freq)))
    #Proportion of NVTs found in NVT-GPSCs or VT-GPSCs
    prop_of_NVTpre <- (pre_counts/(total_NVTpre))*100
    prop_of_NVTpost <- (post_counts/(total_NVTpost))*100
    pre_actual <- sum(subset(pop_tab, Period=="Pre-PCV", select=c(Actual)))
    post_actual <- sum(subset(pop_tab, Period=="Post-PCV13", select=c(Actual)))
    #round to estimated counts
    pop_tab$Actual <- round(pop_tab$Actual)
    #runglm
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
    #Replacement
    extra_cases <- post_inc-pre_inc
    GPSCtype_NVTIRR_pre_PCV13 <- rbind(GPSCtype_NVTIRR_pre_PCV13,c(country, what_type, No_GPSCs, prop_of_NVTpre, prop_of_NVTpost, IRR, confi_lo, confi_up, ps, pre_counts, pre_actual, post_counts, post_actual,pre_inc, post_inc,extra_cases ))
  }
}
colnames(GPSCtype_NVTIRR_pre_PCV13) <- c("Country","GPSC-type", "No. GPSCs","Prop of NVTs pre-PCV","Prop of NVTs pre-PCV","IRR","lower","upper","p", "pre-genomes", "pre-estimated cases","post-genomes","post-estimated cases", "pre-avg-incidence","post-avg-incidence","extra cases")
write.csv(GPSCtype_NVTIRR_pre_PCV13, file ="Tab2_postPCV13_glmIRR_NVT_GPSCtype.csv", row.names = FALSE)
