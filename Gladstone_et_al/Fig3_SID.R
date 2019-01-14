#input data T1-GPSC_dataset tab from Supplementary-T1-T21 in csv format
T8 <- read.csv("T8-GPSC_summaries.csv", header = TRUE, sep =",", stringsAsFactors = FALSE)

#Geographical 1-D plot
geo_dat <- subset(T8, GPSC.type=="Dominant")
pdf("Fig3A_Geographical_simpsons_GPSCs.pdf",width=6,height=6,paper='special')  
boxplot(geo_dat$Country.1.D,geo_dat$Continent.1.D, outline = FALSE, names=c("GPSC (countries)","GPSC (continents)"), ylab="Simpson's 1-D", main = "Geographical diversity of dominant-GPSCs")
stripchart(geo_dat$Country.1.D,  at=1:1, vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
stripchart(geo_dat$Continent.1.D,  at=2:2, vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
dev.off()

#Serotype 1-D pre-PCV plot
Sero_dat <- subset(T8, GPSC.type=="Dominant")
pdf("Fig3B_Serotype_simpsons_GPSCs.pdf",width=6,height=6,paper='special')  
boxplot(Sero_dat$Sero.SDI.1.D.pre.PCV, ylab="Simpson's 1-D", main="Serotype diversity of dominant-GPSCs pre-PCV", xlab="GPSC (serotypes)")
stripchart(Sero_dat$Sero.SDI.1.D.pre.PCV, vertical = TRUE, method = "jitter", add = TRUE, pch = 20)
dev.off()

