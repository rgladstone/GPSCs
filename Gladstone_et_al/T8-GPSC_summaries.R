#input data T1-GPSC_dataset tab from Supplementary-T1-T21 in csv format
SID_data <- read.csv("T1-GPS_dataset.csv", header = TRUE, sep =",", stringsAsFactors = FALSE)
library("vegan")
library("data.table")
library("dplyr")

##Count summaries of GPSC##

#counts of isolates per GPSC
GPSC_counts <- as.data.frame(table(SID_data$GPSC))
colnames(GPSC_counts) <- c("GPSC","Freq")
#counts of isolates collected pre-PCV per GPSC
GPSC_counts_pre <- as.data.frame(table(subset(SID_data, Vaccine_Period=="Pre-PCV",select = c(GPSC))))
colnames(GPSC_counts_pre) <- c("GPSC","Freq.z")
#counts of isolates per GPSC from disease
GPSC_counts_manifest <- as.data.frame(table(subset(SID_data, Clinical_Manifest=="Disease",select = c(GPSC,Clinical_Manifest))))
#counts of isolates collected pre-PCV per GPSC from disease
GPSC_counts_manifest_pre <- as.data.frame(table(subset(SID_data, Clinical_Manifest=="Disease" & Vaccine_Period=="Pre-PCV",select = c(GPSC,Clinical_Manifest))))
#merge together above counts per GPSC into single df
GPSC_counts <- merge(GPSC_counts,GPSC_counts_manifest,by="GPSC", all=TRUE)
GPSC_counts$Clinical_Manifest <- NULL
GPSC_counts <- merge(GPSC_counts,GPSC_counts_pre,by="GPSC", all=TRUE)
GPSC_counts <- merge(GPSC_counts,GPSC_counts_manifest_pre,by="GPSC", all=TRUE)
GPSC_counts$Clinical_Manifest <- NULL
GPSC_counts[is.na(GPSC_counts)] <- 0
colnames(GPSC_counts) <- c("GPSC","Count", "Disease Count", "Pre-PCV Count", "Pre-PCV Disease Count")

# GPSC type (Dominant, Intermediate, Rare)
type <- distinct(subset(SID_data, select = c(GPSC, GPSC.type)))

#subset columns (GPSC,In_Silico_Serotype, Vaccine_Status) for all isolates 
all <- subset(SID_data, select =c(GPSC,In_Silico_Serotype, Vaccine_Status))
#subset columns (GPSC,In_Silico_Serotype, Vaccine_Status) for pre-PCV isolates 
pre <- subset(SID_data, Vaccine_Period=="Pre-PCV", select =c(GPSC,In_Silico_Serotype, Vaccine_Status))

# GPSC by serotype counts table pre-PCV
sero_pre <- table(pre$GPSC,pre$In_Silico_Serotype)
# GPSC by serotype counts table
sero <- table(all$GPSC,all$In_Silico_Serotype)
# GPSC by ST counts table
st <- table(SID_data$GPSC,SID_data$In_Silico_St)
# GPSC by country counts table
countries <- table(SID_data$GPSC,SID_data$Country)
# GPSC by continent counts table
continent <- table(SID_data$GPSC,SID_data$Continent)

# loop through GPSCs 
b <- 0
for (GPSCloop in sort(unique(SID_data$GPSC))){
  #Isolates count for GPSC currently being looped
  GPSC_total <- subset(GPSC_counts, GPSC==GPSCloop, select=c(Count))
  #list serotypes observed in this GPSC
  seros <- unique(subset(all,GPSC==GPSCloop))[,2]
  #If this GPSC was observed pre-PCV then calculate it's pre-PCV serotype diversity using simpson's diversity 1-D index
  if(GPSCloop %in% rownames(sero_pre)){
    seroD1_pre <- diversity(sero_pre[as.character(GPSCloop),], "simpson")
  } else {
    seroD1_pre <- "NA"
  }
  #Calculate the GPSCs country and continental diversity using simpon's diversity 1-D index
  countryD1 <- diversity(countries[as.character(GPSCloop),], "simpson")
  continentD1 <- diversity(continent[as.character(GPSCloop),], "simpson")
  #Count the number of countries and continents observed for the GPSC
  country_count <- length(unique(subset(SID_data, GPSC==GPSCloop, select=c(Country)))[,1])
  continent_count <- length(unique(subset(SID_data, GPSC==GPSCloop, select=c(Continent)))[,1])
  #Isolate counts per county ordered from highest to lowest
  country_df <- as.data.frame(table(subset(SID_data,GPSC==GPSCloop, select =c(Country))))
  country_df  <- country_df [order(-country_df$Freq),,] 
  #country with the most isolates
  country_max <- as.character(country_df[1,1])
  #Isolate count in this GPSC for the country with the most isolates
  max_country_count <- country_df[1,2]
  #Total isolate count for the country with the most isolates
  all_country_count <-sum(SID_data$Country == country_max)
  
  # 2x2 contingency table to test
  country_by2 <- matrix(c(max_country_count, as.numeric(GPSC_total) - max_country_count, all_country_count - max_country_count, nrow(SID_data) - as.numeric(GPSC_total) - (all_country_count - max_country_count)), byrow = TRUE, 2,2)
  colnames(country_by2) <- c("countryX","non-countryX")
  rownames(country_by2) <- c("GPSCX","non-GPSCX")
    # If GPSC has more than 100 isolates (Dominant-GPSCs) test if the number of isolates from the country with the most isolates is greater than expect
  if (GPSC_total>100){
    P_country <- chisq.test(country_by2)$p.value
  } else {  
      P_country <- "NA"  
  }
  
  #Number of isolates per serotype, in this GPSC, ranked
  sero_counts <- table(subset(all,GPSC==GPSCloop, select =c(In_Silico_Serotype)))
  sero_counts_df <- as.data.frame(sero_counts)
  sero_ranked <- sero_counts_df[ order(-sero_counts_df$Freq),,]
  #Number of isolates for the most common serotype and its percentage of that serotype found in this GPSC
  sero_count <- max(sero_counts)
  sero_max <- as.character(sero_ranked[1,1])
  sero_perc <- (sero_count/max(table(subset(all, In_Silico_Serotype==sero_max, select=c(In_Silico_Serotype)))))*100
  
  #The number of isolates expressing serotypes pre-PCV in each of the PCV formulations (licenced and in formulation)
  vaccine_sum_pre <- table(subset(pre,GPSC==GPSCloop, select =c(Vaccine_Status)))
  PCV20 <- (vaccine_sum_pre)[which(names(vaccine_sum_pre)=="PCV20 (Pfizer)")]
  PCV15 <- (vaccine_sum_pre)[which(names(vaccine_sum_pre)=="PCV15 (Merck)")]
  PCV13 <- (vaccine_sum_pre)[which(names(vaccine_sum_pre)=="PCV13")]
  PCV10 <- (vaccine_sum_pre)[which(names(vaccine_sum_pre)=="PCV10")]
  PCV7 <- (vaccine_sum_pre)[which(names(vaccine_sum_pre)=="PCV7")]
  #Percentage of GPSC expressing serotypes included in PCV formulations pre-PCV
  PCV13_preperc <- sum(PCV13,PCV10,PCV7)/sum(vaccine_sum_pre)*100
  if (length(PCV13_preperc)==0 || PCV13_preperc=="NaN"){
    PCV13_preperc <- 0
  }
  PCV10_preperc <- sum(PCV10,PCV7)/sum(vaccine_sum_pre)*100
  if (length(PCV10_preperc)==0 || PCV10_preperc=="NaN"){
    PCV10_preperc <- 0
  } else {
    PCV10_preperc <- sum(PCV10,PCV7)/sum(vaccine_sum_pre)*100
  }
  PCV7_preperc <- PCV7/sum(vaccine_sum_pre)*100
  if (length(PCV7_preperc)==0 || PCV7_preperc=="NaN"){
    PCV7_preperc <- 0
  }
  NVT_preperc <- (sum(vaccine_sum_pre)-sum(PCV7,PCV10,PCV13))/sum(vaccine_sum_pre)*100
  if (length(NVT_preperc)==0 || NVT_preperc=="NaN"){
    NVT_preperc <- 0
  }
  PCV15_preperc <- sum(PCV15,PCV13,PCV10,PCV7)/sum(vaccine_sum_pre)*100
  if (length(PCV15_preperc)==0 || PCV15_preperc=="NaN"){
    PCV15_preperc <- 0
  }
  PCV20_preperc <- sum(PCV20,PCV15,PCV13,PCV10,PCV7)/sum(vaccine_sum_pre)*100
  if (length(PCV20_preperc)==0 || PCV20_preperc=="NaN"){
    PCV20_preperc <- 0
  }
  #The number of isolates expressing serotypes in each of the PCV formulations (licenced and in formulation)
  vaccine_sum <- table(subset(SID_data,GPSC==GPSCloop, select =c(Vaccine_Status))) 
  PCV20 <- (vaccine_sum)[which(names(vaccine_sum)=="PCV20 (Pfizer)")]
  PCV15 <- (vaccine_sum)[which(names(vaccine_sum)=="PCV15 (Merck)")]
  PCV13 <- (vaccine_sum)[which(names(vaccine_sum)=="PCV13")]
  PCV10 <- (vaccine_sum)[which(names(vaccine_sum)=="PCV10")]
  PCV7 <- (vaccine_sum)[which(names(vaccine_sum)=="PCV7")]
  #Percentage of GPSC expressing serotypes included in PCV formulations
  PCV13_perc <- sum(PCV13,PCV10,PCV7)/sum(vaccine_sum)*100
  if (length(PCV13_perc)==0){
    PCV13_perc <- 0
  }
  PCV10_perc <- sum(PCV10,PCV7)/sum(vaccine_sum)*100
  if (length(PCV10_perc)==0){
    PCV10_perc <- 0
  }
  PCV7_perc <- PCV7/sum(vaccine_sum)*100
  if (length(PCV7_perc)==0){
    PCV7_perc <- 0
  }
  NVT_perc <- (sum(vaccine_sum)-sum(PCV7,PCV10,PCV13))/sum(vaccine_sum)*100
  if (length(NVT_perc )== 0){
    NVT_perc <- 0
  }
  PCV15_perc <- sum(PCV15,PCV13,PCV10,PCV7)/sum(vaccine_sum)*100
  if (length(PCV15_perc)==0){
    PCV15_perc <- 0
  }
  PCV20_perc <- sum(PCV20,PCV15,PCV13,PCV10,PCV7)/sum(vaccine_sum)*100
  if (length(PCV20_perc)==0){
    PCV20_perc <- 0
  }
  
  #Number of isolates per ST, in this GPSC
  sts <- unique(subset(SID_data,GPSC==GPSCloop, select =c(In_Silico_St)))[,1]
  st_counts <- table(subset(SID_data,GPSC==GPSCloop, select =c(In_Silico_St)))
  st_counts_df <- as.data.frame(st_counts)
  #most common ST
  st_count <- max(st_counts)
  st_max <- names(st_counts)[which(st_counts==max(st_counts))]
  #mean number of antibiotic classes to which isolates are resistance
  mean_classes <- mean(unique(subset(SID_data,GPSC==GPSCloop, select =c(No_of_classes)))[,1])
  
  #Create output row of summary data for GPSC
  serostD1 <- c(GPSCloop,
                length(seros),seroD1_pre, sero_max, sero_count, sero_perc, paste(as.character(sero_ranked[,1]),collapse=","),PCV7_perc, PCV7_preperc, PCV10_perc, PCV10_preperc,PCV13_perc, PCV13_preperc, PCV15_perc, PCV15_preperc, PCV20_perc, PCV20_preperc, NVT_perc, NVT_preperc, 
                length(sts), st_max[1], st_count,
                continent_count,continentD1,country_count,countryD1,country_max,P_country,mean_classes)
  if (b == 0) {
    results <- serostD1
    b <- 1
  } else {
    results <- rbind(results,serostD1)
  }
}

#remove row names
rownames(results) <- NULL
results <- as.data.frame(results)
colnames(results)[1] <- "GPSC"

#Merge count results with GPSC type
all_results <- merge(GPSC_counts, type, by="GPSC")
#Merge count results with GPSC loop results
all_results <- merge(all_results,results, by="GPSC")
colnames(all_results) <- c("GPSC","Count","Disease Count", "Pre-PCV Count", "Pre-PCV Disease Count", "GPSC type", "Number of Serotypes", "Sero SDI 1-D pre-PCV","Max Serotype", "Max Serotype Count", "% of total serotype X", "Serotype list", "% PCV7", "% PCV7-pre","% PCV10","% PCV10-pre", "% PCV13", "% PCV13-pre","% PCV15","% PCV15-pre","% PCV20","% PCV20-pre", "% PCV13 NVT","% PCV13 NVT-pre", 
                           "Number of STs", "Max ST", "Max ST count","Continent count", "Continent 1-D", "Country count", "Country 1-D", "Country Max", "Chi p max country","Mean_No_resistant_classes")
#Adjust for multiple testing
padj <- subset(all_results, `GPSC type`=="Dominant", c(GPSC, `Chi p max country`)) 
padjd <- p.adjust(as.numeric(as.character(padj$`Chi p max country`)), n= length(padj$`Chi p max country`))
padjd <- cbind(padj,padjd)
padjd$`Chi p max country` <- NULL
colnames(padjd) <- c("GPSC","Chi ajusted p")
#Merge adjustd p-values with all results
final_results <- merge(all_results,padjd, by="GPSC",all=TRUE)
#swap final columns
final_results <- final_results[,c(1:33,35,34)]

#Write file out for T8 of Supplementary tables
fname <- "T8-GPSC_summaries.csv"
fname <- gsub(":","-",fname)                 
fwrite(final_results, file =fname)