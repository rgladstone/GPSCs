library(meta)

#input population size per year from https://github.com/rgladstone/GPSCs/blob/master/lo_et_al/pop_years.csv
pop <- read.csv("pop_years.csv", header = TRUE, sep =",")
pop_post_PCV13 <- subset(pop,Period=="Post-PCV13")
popSA <- subset(pop_post_PCV13,Country=="South Africa", select=c(population,selection))
popIL <- subset(pop_post_PCV13,Country=="Israel", select=c(population,selection))
popUS <- subset(pop_post_PCV13,Country=="USA", select=c(population,selection))
#Average selection term and population post-PCV13
post_selectionSA <- sum(popSA$selection)/nrow(popSA)
post_popSA <- sum(popSA$population)/nrow(popSA)
post_selectionIL <- sum(popIL$selection)/nrow(popIL)
post_popIL <- sum(popIL$population)/nrow(popIL)
post_selectionUS <- sum(popUS$selection)/nrow(popUS)
post_popUS <- sum(popUS)/nrow(popUS)
#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary_v4.csv", header = TRUE, sep =",")

NVT_counts <- droplevels(subset(paper2, Vaccine_Status=="NVT" & Vaccine_Period=="Post-PCV13" & Country %in% c("South Africa", "Israel", "USA"), select = c(GPSC_new, Country)))
NVT_tab <- as.data.frame(unclass(table(NVT_counts$GPSC_new,NVT_counts$Country)))
#estimated cases per country
NVT_tab$SA_estimated <- NVT_tab$`South Africa`/post_selectionSA
NVT_tab$IL_estimated <- NVT_tab$Israel/post_selectionIL
NVT_tab$US_estimated <- NVT_tab$USA/post_selectionUS


GPSC_PCV13_metarate <- matrix(data=NA,nrow=0,ncol=7)
for (gROW in 1:nrow(NVT_tab)){

  casesSA <- NVT_tab[gROW,]$SA_estimated
  casesIL <- NVT_tab[gROW,]$IL_estimated
  casesUS <- NVT_tab[gROW,]$US_estimated
  
  #meta annual incidence rate per 100,000 using log transformation
  meta_IRLN <- metarate(c(casesSA,casesIL,casesUS), c(post_popSA, post_popIL, post_popUS),sm="IRLN")
  #random effect rate
  rmeta_rate <- exp(meta_IRLN$TE.random)*100000
  rmeta_lower <- exp(meta_IRLN$lower.random)*100000
  rmeta_upper <- exp(meta_IRLN$upper.random)*100000
  #fixed effect rate
  fmeta_rate <- exp(meta_IRLN$TE.fixed)*100000
  fmeta_lower <- exp(meta_IRLN$lower.fixed)*100000
  fmeta_upper <- exp(meta_IRLN$upper.fixed)*100000
  GPSC_PCV13_metarate <- rbind(GPSC_PCV13_metarate,c(rownames(NVT_tab)[gROW], rmeta_rate, rmeta_lower, rmeta_upper, fmeta_rate, fmeta_lower, fmeta_upper))
}

colnames(GPSC_PCV13_metarate) <- c("GPSC","RE-metarate","RE-meta lower","RE-meta upper","FE-metarate","FE-meta lower","FE-meta upper")
GPSC_PCV13_metarate <- GPSC_PCV13_metarate[order(as.data.frame(GPSC_PCV13_metarate)$`RE-metarate`,decreasing=TRUE),]
write.csv(GPSC_PCV13_metarate, file ="NVT_pooledmetarate_SA_IL_US.csv", row.names = FALSE)
