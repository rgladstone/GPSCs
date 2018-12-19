#install.packages('VennDiagram')
library(VennDiagram)

#input data T1-GPSC_dataset tab from Supplementary-T1-T20 material in csv format
T1 <- read.csv("T1-GPS_dataset.csv", header = TRUE, sep =",", stringsAsFactors = FALSE)

#Filter down to dominant GPSCs
dom_GPSCs <- subset(T1, GPSC.type=="Dominant" & Vaccine_Period=="Pre-PCV", select=c(GPSC,Vaccine_Status))

PCV_status <- as.data.frame(unclass(table(dom_GPSCs$GPSC, dom_GPSCs$Vaccine_Status)))
#combine NVT and PCV15 and PCV20 VTs as PCV13 NVT
PCV_status$PCV13NVT <- PCV_status$PCV13NVT <- PCV_status[,1]+PCV_status[,4]+PCV_status[,5]
#remove unneccessay columns
PCV_status[c(1,4,5)] <- NULL
#Reorder columns
PCV_status <- PCV_status[c(3,1,2,4)]

#venn areas (number of GPSCs with all combinations of PCV status)
ar1=sum(PCV_status[,1]>0)
ar2=sum(PCV_status[,2]>0)
ar3=sum(PCV_status[,3]>0)
ar4=sum(PCV_status[,4]>0)
ar12=sum(PCV_status[,1]>0 & PCV_status[,2]>0)
ar13=sum(PCV_status[,1]>0 & PCV_status[,3]>0)
ar14=sum(PCV_status[,1]>0 & PCV_status[,4]>0)
ar23=sum(PCV_status[,2]>0 & PCV_status[,3]>0)
ar24=sum(PCV_status[,2]>0 & PCV_status[,4]>0)
ar34=sum(PCV_status[,3]>0 & PCV_status[,4]>0)
ar123=sum(PCV_status[,1]>0 & PCV_status[,2]>0 & PCV_status[,3]>0)
ar124=sum(PCV_status[,1]>0 & PCV_status[,2]>0 & PCV_status[,4]>0)
ar134=sum(PCV_status[,1]>0 & PCV_status[,3]>0 & PCV_status[,4]>0)
ar234=sum(PCV_status[,2]>0 & PCV_status[,3]>0 & PCV_status[,4]>0)
ar1234=sum(PCV_status[,1]>0 & PCV_status[,2]>0 & PCV_status[,3]>0 & PCV_status[,4]>0)

png("Figure4A.png")
licenced_PCVs <- draw.quad.venn(area1=ar1, area2=ar2, area3=ar3,area4=ar4,
                           n12=ar12,n13=ar13,n14=ar14,
                           n23=ar23,n24=ar24,
                           n34=ar34,
                           n123=ar123,n124=ar124,n134=ar134,
                           n234=ar234,
                           n1234=ar1234,
                      category = c("PCV7","PCV10-unique","PCV13-unique","NVTs"), 
                      lty = "blank", fill = c("skyblue", "pink1", "mediumorchid", "orange"))
dev.off()

