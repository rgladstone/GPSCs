#Resistance proportion across dominanat 35 GPSCs
resdata <- read.csv("T1-GPS_dataset.csv", header = TRUE, sep =",")
pre_DM <- subset(resdata, GPSC.type=="Dominant" & Vaccine_Period=="Pre-PCV", select=c(GPSC,Corrected_PEN_MIC_CLSI_meningitis,SXT_predict_by_geno_CLSI,cat,erm_or_mef,tet, MDR,Any_resistant_genotype))
GPSC_counts <- table(pre_DM$GPSC)
tet <-cbind(table(pre_DM$GPSC,pre_DM$tet)[,3],table(pre_DM$GPSC,pre_DM$tet)[,5],table(pre_DM$GPSC,pre_DM$tet)[,6])
tet <- rowSums(tet)
MDR <- as.data.frame(unclass(table(pre_DM$GPSC,pre_DM$MDR)))$yes
penR <- as.data.frame(unclass(table(pre_DM$GPSC,pre_DM$Corrected_PEN_MIC_CLSI_meningitis)))$R
sxtR <- as.data.frame(unclass(table(pre_DM$GPSC,pre_DM$SXT_predict_by_geno_CLSI)))$R
eryR <- as.data.frame(unclass(table(pre_DM$GPSC,pre_DM$erm_or_mef)))$yes
catR <- as.data.frame(unclass(table(pre_DM$GPSC,pre_DM$cat)))$cat1
Results <- cbind(GPSC_counts,
                 penR,
                 sxtR,
                 eryR,
                 catR,
                 tet,
                 MDR,
                 penR/GPSC_counts*100,
                 sxtR/GPSC_counts*100,
                 eryR/GPSC_counts*100,
                 catR/GPSC_counts*100,
                 tet/GPSC_counts*100,
                 MDR/GPSC_counts*100)

colnames(Results) <- c("Total","Pen","SXT","Ery","cat","tet","MDR","Pen.R","Co.tri.R","Ery.R","Chloram.R","Tet.R","MDR.R")

pdf("Fig5_Resistant_proportion_GPSCs.pdf",width=15,height=6,paper='special')  
boxplot(Results[,"Pen.R"],  Results[,"Ery.R"], Results[,"Tet.R"], Results[,"Co.tri.R"], Results[,"Chloram.R"], outline= FALSE, Results[,"MDR.R"], names=c("Penicillin","Erythromycin","Tetracycline","Co-Trimoxazole","Chloramphenicol", "MDR"), main="Antibiotic resistance in dominant-GPSCs pre-PCV", ylab="Percentage resistant per GPSC")
stripchart(Results[,"Pen.R"], at=1:1, vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
stripchart(Results[,"Ery.R"], at=2:2, vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
stripchart(Results[,"Tet.R"], at=3:3,vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
stripchart(Results[,"Co.tri.R"], at=4:4, vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
stripchart(Results[,"Chloram.R"], at=5:5, vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
stripchart(Results[,"MDR.R"], at=6:6, vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
dev.off()
