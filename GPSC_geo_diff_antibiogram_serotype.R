#input data
T1 <- read.csv("/Users/rg9/rg9_documents/GPS/trumps/LID_GPSC/Submission/Github/T1-GPS_dataset.csv", header = TRUE, sep =",")

#Filter down to dominant GPSCs
dom_GPSCs <- subset(T1, GPSC.type=="Dominant" & Vaccine_Period=="Pre-PCV", select=c(GPSC,Country,Antibiogram,In_Silico_Serotype))

#r x c tabluation of antibiogram x GPSC
GPSC_anti_tab <- table(dom_GPSCs$Antibiogram, dom_GPSCs$GPSC)
#test whether the distribution of antibiograms is different across GPSCs, chi-square as too memory intensive for fishers, simulated as small expected values. 
chi_GPSCanti <- chisq.test(GPSC_anti_tab, simulate.p.value = TRUE)
#Result
chi_GPSCanti

#list of dominant GPSCs
GPSCs <- unique(dom_GPSCs$GPSC)
#matrix to add results to
result <- matrix(data=NA,nrow=0,ncol=3)

#loop through each lineage and test (chi.sq) whether the distribution of A) Serotypes and B) Antibiograms is different across countries.
for (lineage in GPSCs){
  #subset on GPSC
  GPSC_subset <- subset(dom_GPSCs, GPSC==lineage, select = c(Country, In_Silico_Serotype, Antibiogram))
  #make factors/drop unused factors
  GPSC_subset$Country <- factor(GPSC_subset$Country)
  GPSC_subset$Antibiogram <- factor(GPSC_subset$Antibiogram)
  GPSC_subset$In_Silico_Serotype <- factor(GPSC_subset$In_Silico_Serotype)
  #r x c tabluation of antibiogram x country
  anti_tab <- table(GPSC_subset$Antibiogram, GPSC_subset$Country)
  #chi-square as too memory intensive for fishers, simulated as small expected values.
  chi_anti <- chisq.test(anti_tab, simulate.p.value = TRUE)
  #r x c tabluation of serotype x country
  sero_tab <- table(GPSC_subset$In_Silico_Serotype, GPSC_subset$Country)
  #chi-square as too memory intensive for fishers, simulated as small expected values.
  chi_sero <- chisq.test(sero_tab, simulate.p.value = TRUE)
  #results
  result <- rbind(result,c(lineage,chi_sero$p.value,chi_anti$p.value))
}

#Adjust for multiple testing

result <- cbind(result, p.adjust(as.numeric(result[,2]), method = "BH", n=length(result[,2])), p.adjust(as.numeric(result[,3]), method = "BH", n=length(result[,3])))
colnames(result) <- c("GPSC","Sero_p", "Anti_p", "Sero_adjp", "Anti_adjp")
#Count how many GPSCs had p<0.01 strict cut off due to small values which will make estimates poor 
results <- cbind(length(which(as.numeric(result[,4]) < 0.01)), length(which(as.numeric(result[,5]) < 0.01)))
colnames(results) <- c("No. dominant GPSCs with Sig dif distribution of serotype across countries", "No. dominant GPSCs with Sig dif distribution of antibiogram across countries")
rownames(results) <- "X/35"
#print final result
results
