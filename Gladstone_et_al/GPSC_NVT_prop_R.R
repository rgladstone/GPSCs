#input data T1-GPSC_dataset tab from Supplementary-T1-T21 in csv format
res_data <- read.csv("T1-GPS_dataset.csv", header = TRUE, sep =",", stringsAsFactors = FALSE)

#PCV VT proportion pre vs post for GPSCs that contain PCV7-VT and PCV13-NVTs
Results <- matrix(data=NA,nrow=0,ncol=8)

dominant_GPSC <-  unique(subset(res_data, GPSC.type=="Dominant", select =c(GPSC)))

for (x in dominant_GPSC[,1]){
  by2 <- as.data.frame(unclass(table(subset(res_data, GPSC==x, select =c(Corrected_PEN_MIC_CLSI_meningitis, PCV13_Status)))))
  if (length(by2[1,])>1 & length(by2[,1])>1 ){
    if(sum(by2[1,1],by2[2,1])>=10 & sum(by2[1,2],by2[2,2])>=10){
    prop <- fisher.test(by2)$p.value
    } else {
    prop <- "NA"
    }
    Results <- rbind(Results,c(x, by2[1,1],by2[2,1], by2[1,1]/sum(by2[1,1],by2[2,1]), by2[1,2], by2[2,2], by2[1,2]/sum(by2[1,2],by2[2,2]), prop))
  }
}
colnames(Results) <- c("GPSC","NVT R", "NVT S", "NVT %R", "VT R", "VT S", "VT %R", "p-value")
Results <- cbind(Results, p.adjust(Results[,8], method = "BH"))
colnames(Results)[9] <- "adjusted p-value"
write.csv(Results, file ="T14-NVT_penicillin_resistance.csv", row.names = FALSE)
