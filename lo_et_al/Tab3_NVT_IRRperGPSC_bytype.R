require(tidyverse)

#input population size per year from https://github.com/rgladstone/GPSCs/blob/master/lo_et_al/pop_years.csv
pop <- read.csv("pop_years.csv", header = TRUE, sep =",")
#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary_v4.csv", header = TRUE, sep =",")

#######Calculate NVT IRR per GPSC############
GPSC_NVTIRR_pre_PCV13 <- matrix(data=NA,nrow=0,ncol=9)
GPSC_NVTIRR_pre_PCV7 <- matrix(data=NA,nrow=0,ncol=9)
GPSC_NVTIRR_PCV7_PCV13 <- matrix(data=NA,nrow=0,ncol=7)

for (country in unique(pop$Country)){
  #Determine GPSC type per country based on VT proportion (<50% or >=50%) in first period observed
  GPSC_props <- matrix(data=NA,nrow=0,ncol=3)
  colnames(GPSC_props) <- c("GPSC","Proportion VT", "type")
  coun_p2 <- subset(paper2, Country==country, select=c(GPSC,Vaccine_Period,Vaccine_Status,In_Silico_Serotype,Country,Year))
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
                    Country==country & Vaccine_Period=="Pre-PCV"& Vaccine_Status=="NVT",
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
      #adjust genome count by case selection
      pop_tab$Actual <- pop_tab$Freq/pop_tab$selection
      #round to estimated counts
      pop_tab$Actual <- round(pop_tab$Actual)
      #Only test GPSCs with >10 estimated NVT cases
      if  (sum(pop_tab$Actual)>10){
        if (sum(subset(pop_tab, Period==post)['Actual'])==0){
          pop_tab$Actual <- pop_tab$Actual+1
        }
        if (sum(subset(pop_tab, Period=="Pre-PCV")['Actual'])==0){
          pop_tab$Actual <- pop_tab$Actual+1
        }
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
        if (post == "Post-PCV7"){
          GPSC_NVTIRR_pre_PCV7 <- rbind(GPSC_NVTIRR_pre_PCV7,c(country, cluster, what_type, sum(subset(pop_tab, Period=="Pre-PCV")['Actual']), sum(subset(pop_tab, Period==post)['Actual']), IRR, confi_lo, confi_up, ps))
        } else {
          GPSC_NVTIRR_pre_PCV13 <- rbind(GPSC_NVTIRR_pre_PCV13,c(country, cluster, what_type, sum(subset(pop_tab, Period=="Pre-PCV")['Actual']), sum(subset(pop_tab, Period==post)['Actual']), IRR, confi_lo, confi_up, ps))
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
    pop_tab$Actual <- pop_tab$Freq/pop_tab$selection
    pop_tab$Actual <- round(pop_tab$Actual)
    if  (sum(pop_tab$Actual)>10){
      if (sum(subset(pop_tab, Period==post)['Actual'])==0){
        pop_tab$Actual <- pop_tab$Actual+1
      }
      if (sum(subset(pop_tab, Period=="Post-PCV7")['Actual'])==0){
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
      GPSC_NVTIRR_PCV7_PCV13 <- rbind(GPSC_NVTIRR_PCV7_PCV13,c(country, cluster, what_type, IRR, confi_lo, confi_up, ps))
    }
  }
}



GPSC_NVTIRR_pre_PCV7 <- cbind(GPSC_NVTIRR_pre_PCV7,p.adjust(GPSC_NVTIRR_pre_PCV7[,9], method = "BH", n=length(GPSC_NVTIRR_pre_PCV7[,9])))
GPSC_NVTIRR_pre_PCV13 <- cbind(GPSC_NVTIRR_pre_PCV13,p.adjust(GPSC_NVTIRR_pre_PCV13[,9], method = "BH", n=length(GPSC_NVTIRR_pre_PCV13[,9])))
GPSC_NVTIRR_PCV7_PCV13 <- cbind(GPSC_NVTIRR_PCV7_PCV13,p.adjust(GPSC_NVTIRR_PCV7_PCV13[,7], method = "BH", n=length(GPSC_NVTIRR_PCV7_PCV13[,7])))
colnames(GPSC_NVTIRR_pre_PCV7) <- c("Country","GPSC","GPSC_type","Pre-PCV cases", "Post cases", "IRR","lower","upper","p", "adj.p")
colnames(GPSC_NVTIRR_pre_PCV13) <- c("Country","GPSC","GPSC_type","Pre-PCV cases", "Post cases","IRR","lower","upper","p", "adj.p")
colnames(GPSC_NVTIRR_PCV7_PCV13) <- c("Country","GPSC","GPSC_type", "IRR","lower","upper","p", "adj.p")
write.csv(GPSC_NVTIRR_pre_PCV7, file ="postPCV7_glmIRR_NVT_pseudo_adjfreq.csv", row.names = FALSE)
write.csv(GPSC_NVTIRR_pre_PCV13 , file ="postPCV13_glmIRR_NVT_pseudo_adjfreq.csv", row.names = FALSE)
write.csv(GPSC_NVTIRR_PCV7_PCV13 , file ="postPCV7_postPCV13_glmIRR_NVT_pseudo_adjfreq.csv", row.names = FALSE)

pre_v_PCV7_pre_v_PCV13 <- merge(GPSC_NVTIRR_pre_PCV7[,c(1:3,6,10)],GPSC_NVTIRR_pre_PCV13[,c(1,2,6,10)], by=c("Country","GPSC"), all=TRUE)
pre_v_PCV7_pre_v_PCV13_PCV7_v_PCV13 <- merge(pre_v_PCV7_pre_v_PCV13,GPSC_NVTIRR_PCV7_PCV13[,c(1,2,4,8)], by=c("Country","GPSC"), all=TRUE)
colnames(pre_v_PCV7_pre_v_PCV13_PCV7_v_PCV13) <- c("Country","GPSC","GPSC_type","IRR-pre_v_PCV7","adj.p","IRR-pre_v_PCV13","adj.p","IRR-PCV7_v_PCV13","adj.p")
write.csv(pre_v_PCV7_pre_v_PCV13_PCV7_v_PCV13, file ="Tab3_multiperiod_summary_glmIRR_NVT_pseudo_adjfreq.csv", row.names = FALSE)

