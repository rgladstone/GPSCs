library(RColorBrewer)

#input data T1-GPSC_dataset tab from Supplementary-T1-T20 material in csv format
T1 <- read.csv("T1-GPS_dataset.csv", header = TRUE, sep =",", stringsAsFactors = FALSE)

#Ranked counts of GPSCs
x <- as.vector(sort(as.data.frame(unclass(table(T1$GPSC)))[,1], decreasing=TRUE))

#Categorise GPSCs Dominant n>100, Intermediate n>10<=100, Rare n<10
GPSCcat <- cut(x,c(800,100,10,0))
#Set category colours
cols <- brewer.pal(n=3, name="Set1")
cols_T1 <- cols[GPSCcat]
#Write Figure S4
png("FigureS4_GPSC_types.png")
plot(x, xlab="GPSCs", xaxt='n', ylab="Isolate count", main="GPSCs ranked by representation in GPS (n=13,454)", col=cols_T1)
legend(416,800,legend = unique(GPSC_dat$Classification), col = cols,pch=1)
axis(side=2,at=seq(0, 800, by = 100)) 
abline(v=which.max(elbow(x)))

