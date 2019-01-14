#input data T1-GPSC_dataset tab from Supplementary-T1-T21 in csv format
T8 <- read.csv("T8-GPSC_summaries.csv", header = TRUE, sep =",", stringsAsFactors = FALSE)

GPSC_data <- subset(T8, GPSC.type=="Dominant")
pdf("Fig4B_perc_PCV13VT_perGPSC.pdf",width=6,height=6,paper='special')  
boxplot(GPSC_data$X..PCV13.pre, ylab="Percentage PCV13-VT", main="PCV13-VT percentage pre-PCV", xlab="Dominant-GPSC")
stripchart(GPSC_data$X..PCV13.pre, vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
dev.off()
