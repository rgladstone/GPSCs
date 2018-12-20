library(ggplot2)

#input data T1-GPSC_dataset tab from Supplementary-T1-T20 material in csv format
T1 <- read.csv("T1-GPS_dataset.csv", header = TRUE, sep =",", stringsAsFactors = FALSE)

#list of dominant-GPSCs
dom <- unique(subset(T1, GPSC.type=="Dominant", select = c(GPSC)))[,1]
#Create results file
dists <- data.frame()
#loop through GPSCs
for (cluster in dom){
  pairSNP_cluster <- read.csv(paste("GPSC", cluster, "_noOGref.recombremoved_pairsnp.tsv", sep=""), header = FALSE, sep ="\t")
  #lanes as rownames then conver to matrix
  rownames(pairSNP_cluster) <- pairSNP_cluster$V1
  pairSNP_cluster$V1 <- NULL
  pairSNP_cluster <- as.matrix(pairSNP_cluster)
  #remove self distances
  diag(pairSNP_cluster) <- NA
  #remove duplicated distances in lower triangle
  pairSNP_cluster[lower.tri(pairSNP_cluster)] <- NA
  #collapse matrix into vector of all distances
  cluster_distances <- as.vector(pairSNP_cluster)
  #drop NA values (self, lower triangle)
  cluster_distances <- cluster_distances[!is.na(cluster_distances)]
  #convert to df
  cluster_distances <- as.data.frame(cluster_distances)
  #Add column with the GPSC name per isolate
  cluster_distances$GPSC <- cluster
  #add rows for GPSC to results file
  dists <- rbind(dists, cluster_distances)
}

colnames(dists) <- c("pairSNPs","GPSC")
dists$GPSC <- factor(dists$GPSC)
#plot Figure S5
ggplot(data = dists, aes(x=GPSC, y=pairSNPs)) + geom_violin(size=0.1,fill='#FF5733') + ggtitle("Pairwise SNP distances within GPSCs after recombination removed") + ylab("mutations")
ggsave("/Users/rg9/rg9_documents/GPS/trumps/LID_GPSC/Paper1 version4/recombremoved_SNP_distances_perGPSC.png", plot = last_plot(), device = "png", path = NULL, scale = 1, width = 60, height = 15, units = "cm", dpi = 300, limitsize = TRUE)
