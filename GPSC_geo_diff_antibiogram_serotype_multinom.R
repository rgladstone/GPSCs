#input data T1-GPSC_dataset tab from Supplementary-T1-T20 material in csv format
T1 <- read.csv("/Users/rg9/rg9_documents/GPS/trumps/LID_GPSC/Submission/Github/T1-GPS_dataset.csv", header = TRUE, sep =",")

pre_dom <- subset(T1, GPSC.type=="Dominant" & Vaccine_Period=="Pre-PCV", select=c(Antibiogram))

#top antibiograms with >15 in pre-PCV dom-GPSCs, sample retention n=4037/4122
top_anti <- rownames(subset(as.data.frame(unclass(table(pre_dom$Antibiogram))),  as.data.frame(unclass(table(pre_dom$Antibiogram)))[,1] >=15))
#Filter down to dominant GPSCs
dom_GPSCs <- subset(T1, GPSC.type=="Dominant" & Vaccine_Period=="Pre-PCV" & Antibiogram %in% top_anti, select=c(GPSC,Country,Antibiogram,In_Silico_Serotype))

#Convert to factors
dom_GPSCs$Antibiogram <- as.factor(dom_GPSCs$Antibiogram)
dom_GPSCs$GPSC <- as.factor(dom_GPSCs$GPSC)

#set reference to fully susceptable antibiogram
dom_GPSCs$Antibiogram <- relevel(dom_GPSCs$Antibiogram, ref="1")

#multinomial regression for GPSC as a predictor of antibiogram
anti_mnom <- multinom(Antibiogram ~ GPSC, data=dom_GPSCs, maxit=1000)
#converged?
conv <- anti_mnom$convergence == 0
#Predict
cm <- table(predict(anti_mnom),dom_GPSCs$Antibiogram)
#Misclass error
m_err <- 1-sum(diag(cm))/sum(cm)
# use sum coding, necessary to make type III LR tests valid
set_sum_contrasts()
#likelihood ratio test
anova_lrt <- Anova(anti_mnom,type="III")
#better than null?
p <- anova_lrt$`Pr(>Chi)`[2]

#GPSC as predictor for serotype
sero_dom <- subset(T1, GPSC.type=="Dominant" & Vaccine_Period=="Pre-PCV", select=c(In_Silico_Serotype))

#top serotypes with >15 in pre-PCV dom-GPSCs, sample retention n=/4122
top_sero <- rownames(subset(as.data.frame(unclass(table(pre_dom$In_Silico_Serotype))),  as.data.frame(unclass(table(pre_dom$In_Silico_Serotype)))[,1] >=15))
#Filter down to dominant GPSCs
dom_GPSCs <- subset(T1, GPSC.type=="Dominant" & Vaccine_Period=="Pre-PCV" & In_Silico_Serotype %in% top_sero, select=c(GPSC,Country,Antibiogram,In_Silico_Serotype))

#Convert to factors
dom_GPSCs$In_Silico_Serotype <- droplevels(dom_GPSCs$In_Silico_Serotype)
dom_GPSCs$GPSC <- as.factor(dom_GPSCs$GPSC)

#multinomial regression for GPSC as a predictor of serotype
sero_mnom <- multinom(In_Silico_Serotype ~ GPSC, data=dom_GPSCs, maxit=1000)
#converged?
conv <- sero_mnom$convergence == 0
#Predict
cm_sero <- table(predict(sero_mnom),dom_GPSCs$In_Silico_Serotype)
#Misclass error
m_err_sero <- 1-sum(diag(cm_sero))/sum(cm_sero)
# use sum coding, necessary to make type III LR tests valid
set_sum_contrasts()
#likelihood ratio test
anova_lrt_sero <- Anova(sero_mnom,type="III")
#better than null?
p_sero <- anova_lrt_sero$`Pr(>Chisq)`

#reset to dominant GPSCs pre-PCV
dom_GPSCs <- subset(T1, GPSC.type=="Dominant" & Vaccine_Period=="Pre-PCV", select=c(GPSC,Country,Antibiogram,In_Silico_Serotype))

#list of dominant GPSCs
GPSCs <- unique(dom_GPSCs$GPSC)
#matrix to add results to
result <- matrix(data=NA,nrow=0,ncol=7)

#loop through each lineage and determine whether country is a better than null predictor of A) Antibiogram B) serotype and the proportion of variance it explains
for (lineage in GPSCs){
  GPSC_subset <- subset(dom_GPSCs, GPSC==lineage, select = c(Country, In_Silico_Serotype, Antibiogram))
  #make factors/drop unused factors
  GPSC_subset$Country <- factor(GPSC_subset$Country)
  GPSC_subset$Antibiogram <- factor(GPSC_subset$Antibiogram)
  GPSC_subset$In_Silico_Serotype <- factor(GPSC_subset$In_Silico_Serotype)
  #set reference to susceptable (if present), or to the most common profile
  if (1 %in% levels(GPSC_subset$Antibiogram)){
    GPSC_subset$Antibiogram <- relevel(GPSC_subset$Antibiogram, ref=1)
  } else {
    GPSC_subset$Antibiogram <- relevel(GPSC_subset$Antibiogram, ref=levels(GPSC_subset$Antibiogram)[1]) }                                   
  
  #country as a predictor of antibiogram in multinomial
  anti_country_mnom <- multinom(Antibiogram ~ Country, data=GPSC_subset, maxit=1000)
  #converged?
  convanti <- anti_country_mnom$convergence == 0
  #predict
  country_cm <- table(predict(anti_country_mnom),GPSC_subset$Antibiogram)
  #Misclass error
  country_m_err <- 1-sum(diag(country_cm))/sum(country_cm)
  # use sum coding, necessary to make type III LR tests valid
  set_sum_contrasts()
  #likelihood ratio test
  country_anova_lrt <- Anova(anti_country_mnom,type="III")
  #better than null?
  loop_p <- country_anova_lrt$`Pr(>Chisq)`
  
  #Country as a predictor of serotype in multinomial
  if (length(levels(GPSC_subset$In_Silico_Serotype)) >1){
    sero_country_mnom <- multinom(In_Silico_Serotype ~ Country, data=GPSC_subset, maxit=1000)
    #converged?
    convsero <- sero_country_mnom$convergence == 0
    #confusion matrix
    sero_cm <- table(predict(sero_country_mnom),GPSC_subset$In_Silico_Serotype)
    #Misclass error
    sero_m_err <- 1-sum(diag(sero_cm))/sum(sero_cm)
    # use sum coding, necessary to make type III LR tests valid
    set_sum_contrasts()
    #likelihood ratio test
    sero_anova_lrt <- Anova(sero_country_mnom,type="III")
    #better than null?
    loop_sero_p <- sero_anova_lrt$`Pr(>Chisq)`
  } else {
    sero_m_err <- "NA"
    loop_sero_p <- "NA"
    convsero <- "NA"
  }
  
  #results
  result <- rbind(result,c(lineage, as.character(convanti), country_m_err, loop_p, as.character(convsero), sero_m_err, loop_sero_p))
}

#Adjust for multiple testing
result <- cbind(result, p.adjust(as.numeric(result[,4]), method = "BH"),p.adjust(as.numeric(result[,7]), method = "BH"))
colnames(result) <- c("GPSC", "countryanti_convergance", "countryanti_misclass", "countryanti_multinom_p", "countrysero_convergance","countrysero_misclass", "countrysero_multinom_p", "countryanti_multinom_adjp", "countrysero_multinom_adjp")
#Count dominant-GPSCs where country was sig pred of antibiogram or serotype and report their mean and stdev of _misclass
antisub <- subset(as.data.frame(result), as.numeric(result[,8]) < 0.05, select = c(countryanti_misclass))[,1]
serosub <- subset(as.data.frame(result), as.numeric(result[,9]) < 0.05, select = c(countrysero_misclass))[,1]
results <- cbind(length(which(as.numeric(result[,8]) < 0.05)),
                 mean(as.numeric(as.matrix(antisub))),
                 sd(as.numeric(as.matrix(antisub))),
                 length(which(as.numeric(result[,9]) < 0.05)),
                 mean(as.numeric(as.matrix(serosub))),
                 sd(as.numeric(as.matrix(serosub))))
colnames(results) <- c("N/35 dom-GPSCs country sigpred for antibiogram", "mean misclass country sigpred for antibiogram", "stdev misclass country sigpred for antibiogram",
                       "N/35 dom-GPSCs country sigpred for serotype","mean misclass country sigpred for serotype", "stdev misclass country sigpred for serotype")
#print final result
results
