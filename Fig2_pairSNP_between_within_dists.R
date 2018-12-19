library(reshape2)
library(ggplot2)
#input data T1-GPSC_dataset tab from Supplementary-T1-T20 material in csv format
T1 <- read.csv("T1-GPS_dataset.csv", header = TRUE, sep =",", stringsAsFactors = FALSE)

#SNP distance matrix from https://github.com/gtonkinhill/pairsnp run on Roary core gene alignment reduced to SNP sites
pairSNPs <- read.csv("core_snp_dists.tsv", header = FALSE, sep ="\t")
pairSNP_mat <- as.matrix(pairSNPs[,-1])
rownames(pairSNP_mat) <- pairSNPs[,1]
colnames(pairSNP_mat) <- rownames(pairSNP_mat)

#lanes and correspondings GPSCs
GPSC <- subset(T1, select=c(ID,GPSC))
#lanes in order found in pairSNPs
lanes <- as.data.frame(pairSNPs[,1])
colnames(lanes) <- "ID"
#map clusters back onto lanes IDs in the order found in pairSNPs
clusters <- merge(lanes,GPSC, by.x = "ID", all.x = TRUE)
clusters_mat <- as.matrix(clusters[,-1])
rownames(clusters_mat) <- clusters[,1]
#set n=8 isolates later excluded due to bad data with cluster NA to cluster "0"
clusters_mat[is.na(clusters_mat)] <- 0 
#lane and cluster data from matrix into df format
cluster_df <- as.data.frame(clusters_mat[,1])
rownames(cluster_df) <- NULL
colnames(cluster_df) <- c("cluster")
#List of GPSCs to loop through
GPSCs <- sort(as.numeric(unique(clusters_mat[,1])))

##seperate between and within distances##
#create copy of matrix
pairSNP_mat2 <- cbind(pairSNP_mat)
#Remove self distances
diag(pairSNP_mat2) <- NA
#Remove lower triangular of matrix (duplicated distances)
pairSNP_mat2[lower.tri(pairSNP_mat2)] <- NA

#Create within distances results file
within_dists <- vector("list",length = 520)

#Collect distances within each GPSC and replace distances with NA in matrix once collected
for (i in GPSCs){
  sub_cluster <- subset(cluster_df, cluster==i)
  if (nrow(sub_cluster) >1){
    pairs <- t(combn(rownames(sub_cluster),2))
    pairs <- apply(pairs,2,as.numeric)
    #exclude distance isolates in cluster 0 (bad data), before calculating distances
    if (i ==0){
      pairSNP_mat2[pairs] <- NA
      pairSNP_mat2[as.numeric(rownames(sub_cluster)),] <- NA
      pairSNP_mat2[,as.numeric(rownames(sub_cluster))] <- NA
    } else{
      if (is.null(dim(pairs))){
        pairs <-  t(as.matrix(pairs))
        within_dists[[match(i,GPSCs)]] <- pairSNP_mat2[pairs]
        pairSNP_mat2[pairs] <- NA
      } else{
        within_dists[[match(i,GPSCs)]] <- pairSNP_mat2[pairs]
        pairSNP_mat2[pairs] <- NA
      }
    }
  }
}

#vector of all distances from within GPSCs
all_within_dists <- unlist(within_dists, recursive = TRUE)

#remove NA distances (only one isolate in GPSC no distances recorded)
all_within_dists <- all_within_dists[!is.na(all_within_dists)]

#vector of all remaining distances representing between distances
between_values <- melt(pairSNP_mat2)
#remove NA distances (within GPSC dists, excluded isolates distances, self distances, lower triangle)
between_dists <- between_values$value[!is.na( between_values$value)]

#Within distances in df
within_df <- as.data.frame(all_within_dists)
within_df$Type <- "Within"
colnames(within_df) <- c("SNPdist","Type")
#Between distances in df
between_df <- as.data.frame(between_dists)
between_df$Type <- "Between"
colnames(between_df) <- c("SNPdist","Type")

#within and between in one df
SNP_dists <- rbind(within_df,between_df)

#Plot Figure 2
ggplot(data = as.data.frame(SNP_dists), aes(x=Type, y=SNPdist)) + geom_violin() + ggtitle("Pairwise core SNP distances between and within GPSCs") + ylab("SNP distances")
#Save Figure 2
ggsave("/Users/rg9/rg9_documents/GPS/trumps/LID_GPSC/Paper1 version4/coreSNP_distances_between_vs_within.png", plot = last_plot(), device = "png", path = NULL, scale = 1, width = 20, height = 20, units = "cm", dpi = 300, limitsize = TRUE)
