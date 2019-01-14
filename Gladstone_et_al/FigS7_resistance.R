library(ggplot2)
library(stringr)

#input data T1-GPSC_dataset tab from Supplementary-T1-T21 in csv format
T1 <- read.csv("T1-GPS_dataset.csv", header = TRUE, sep =",", stringsAsFactors = FALSE)

resist_cols <- subset(T1, select=c(Corrected_CLSI_Pen.SXT.cat.ery.tet, MDR))
resist_cols$MDR[resist_cols$MDR == "no"] <- "S"
resist_cols$MDR[resist_cols$MDR == "yes"] <- "R"
resist_cols$Any[grepl("R", resist_cols$Corrected_CLSI_Pen.SXT.cat.ery.tet) == TRUE] <- "R"
resist_cols$Any[!grepl("R", resist_cols$Corrected_CLSI_Pen.SXT.cat.ery.tet) == TRUE & grepl("I", resist_cols$Corrected_CLSI_Pen.SXT.cat.ery.tet) == TRUE ] <- "I"
resist_cols$Any[is.na(resist_cols$Any)] <- "S"
resist <- str_split_fixed(resist_cols$Corrected_CLSI_Pen.SXT.cat.ery.tet, "_", 5)
resist <- as.data.frame(cbind(resist,resist_cols$MDR, resist_cols$Any))
colnames(resist) <- c("Penicillin","Co-trimoxazole", "Chloramphenicol","Erythromycin","Tetracyline", "MDR", "Any")
results <-matrix(nrow = 0, ncol=2)
colnames(results) <- c("Antibiotic","Status")
for (rescol in colnames(resist)){
  resres <- cbind(rescol,as.data.frame(resist[rescol]))
  colnames(resres) <- c("Antibiotic", "Status")
  results <- rbind(results, resres)
}

ggplot(results, aes(Antibiotic, fill = factor(Status, levels=c("S","I","R")))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(y = "Percentage") + 
  scale_fill_manual(breaks=c("S","I","R"),labels=c("Sensitive","Intermediate","Resistant"),name = "Status",values = c("#2dc937", "#e7b416", "#cc3232"))
ggsave("FigS7_resistance_prevalence.png", plot = last_plot(), device = "png", path = NULL, scale = 1, width = 25, height = 15, units = "cm", dpi = 300, limitsize = TRUE)

