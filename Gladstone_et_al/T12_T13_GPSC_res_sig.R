library("data.table")

#input data T1-GPSC_dataset tab from Supplementary-T1-T21 in csv format
T1 <- read.csv("T1-GPS_dataset.csv", header = TRUE, sep =",")

#GPSCs fishers
 
pre_pen_collection <-  table(subset(T1, Vaccine_Period=="Pre-PCV" , c(Corrected_PEN_MIC_CLSI_meningitis)))

pre_pen_GPSC <-  table(subset(T1, Vaccine_Period=="Pre-PCV" & GPSC.type=="Dominant", c(GPSC,Corrected_PEN_MIC_CLSI_meningitis)))
 
pen_GPSCresults <- matrix(data=NA,nrow=0,ncol=4)
 
 for (x in row.names(pre_pen_GPSC)){
   fishers <- "NA"
   pen_GPSC_by2 <- matrix(c(pre_pen_GPSC[x,2],pre_pen_GPSC[x,3],pre_pen_collection[2]-pre_pen_GPSC[x,2],pre_pen_collection[3]-pre_pen_GPSC[x,3]), byrow = TRUE, 2,2)
   colnames(pen_GPSC_by2) <- c("R","S")
   rownames(pen_GPSC_by2) <- c("GPSC","non-GPSC")
   pen_GPSC_by2 <- as.table(pen_GPSC_by2)
   if (sum(pen_GPSC_by2[1,])>=10){
     fishers <- fisher.test(pen_GPSC_by2, alternative = "greater")
     pen_GPSCresults <- rbind(pen_GPSCresults,c(x,pen_GPSC_by2[1,1],pen_GPSC_by2[1,2],as.numeric(fishers$p.value)))
   } else {
     pen_GPSCresults <- rbind(pen_GPSCresults,c(x,pen_GPSC_by2[1,1],pen_GPSC_by2[1,2],"NA"))
   }   
 }
 colnames(pen_GPSCresults) <- c("GPSC","R","S","p-value")

 pen_GPSCresults <- as.data.frame(pen_GPSCresults,stringsAsFactors = FALSE)
 pen_GPSCresults$adjp <- p.adjust(as.numeric(pen_GPSCresults$`p-value`), method="BH")
 fwrite(pen_GPSCresults, file ="T12-GPSC_penicillin.csv")
 
 pre_MDR_collection <-  table(subset(T1, Vaccine_Period=="Pre-PCV", c(MDR)))
 
 pre_MDR_GPSC <-  table(subset(T1, Vaccine_Period=="Pre-PCV"& GPSC.type=="Dominant" , c(GPSC,MDR)))
 
 MDR_GPSCresults <- matrix(data=NA,nrow=0,ncol=4)
 
 for (x in row.names(pre_MDR_GPSC)){
   fishers <- "NA"
   MDR_GPSC_by2 <- matrix(c(pre_MDR_GPSC[x,3],pre_MDR_GPSC[x,2],pre_MDR_collection[3]-pre_MDR_GPSC[x,3],pre_MDR_collection[2]-pre_MDR_GPSC[x,2]), byrow = TRUE, 2,2)
   colnames(MDR_GPSC_by2) <- c("Yes","No")
   rownames(MDR_GPSC_by2) <- c("GPSC","non-GPSC")
   MDR_GPSC_by2 <- as.table(MDR_GPSC_by2)
   if (sum(MDR_GPSC_by2[1,])>=10){
     fishers <- fisher.test(MDR_GPSC_by2, alternative = "greater")
     MDR_GPSCresults <- rbind(MDR_GPSCresults, c(x,MDR_GPSC_by2[1,1],MDR_GPSC_by2[1,2],as.numeric(fishers$p.value)))
   } else {
     MDR_GPSCresults <- rbind(MDR_GPSCresults, c(x,MDR_GPSC_by2[1,1],MDR_GPSC_by2[1,2],"NA"))
   }   
 }
 colnames(MDR_GPSCresults) <- c("GPSC","Y","N","p-value")
 
 MDR_GPSCresults <- as.data.frame(MDR_GPSCresults,stringsAsFactors = FALSE)
 colnames(MDR_GPSCresults) <- c("GPSC","Yes","No","p-value")
 MDR_GPSCresults$adjp <- p.adjust(as.numeric(MDR_GPSCresults$`p-value`), method="BH")
 fwrite(MDR_GPSCresults, file ="T13-GPSC_MDR.csv")
 