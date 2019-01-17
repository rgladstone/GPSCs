library("vegan")

#input data T2-GPSC assignment dataset tab from Supplementary-T1-T21 in csv format
gpsc_accum <- read.csv("T2-GPSC_assignment_dataset.csv", header = TRUE, sep =",")


##Species richness##

#Format data for species richness analysis
spec_rich <- as.matrix(table(gpsc_accum$Taxon,gpsc_accum$GPSC))
#replace positive frequencies with 1
spec_rich <- 1*(spec_rich>0)
#Estimate species richness
est <- specpool(spec_rich)
print(est)

##GPSC sampling and accumulation##

#Create matrix for storing results
random_input <- data.frame(sort(unique(gpsc_accum$GPSC)))
row.names(random_input) <- random_input$sort.unique.gpsc_accum.GPSC..
random_input[1] <- NULL

#Create community dataset input for specaccum()
#Select a random 380 samples from each of the 12 countries where n>380, merge country samples on GPSCs rows to store GPSC counts per country
for (country in c("South Africa", "BRAZIL","USA", "UK", "China", 
                  "Peru", "Israel", "Malawi", "Nepal", "Netherlands","Thailand", "The Gambia" )){
  country_set <- subset(gpsc_accum,Country==country, select = c(Taxon,GPSC))
  random_380 <-  country_set[sample(nrow(country_set),380),]
  format_380 <- as.data.frame(table(random_380$GPSC))
  row.names(format_380) <- format_380$Var1
  format_380[1] <- NULL
  colnames(format_380) <- paste(country, "freq", sep=".")
  random_input <- merge(random_input,format_380,by="row.names",all=TRUE)
  rownames(random_input) <- random_input$Row.names
  random_input$Row.names <- NULL
}
rowSums(!is.na(random_input))
#replace NA with 0 and transpose
random_input[is.na(random_input)] <- 0
no_GPSCs <- table(rowSums(random_input)>0)
random_input <- t(random_input)

#run 100 random permutations of the community dataset, where sites (countries) are added in random order to determine GPSC accumulation curve
accum_country_380 <- specaccum(random_input, method="random", permutations = 100)
png(file = "GPSC_accum_12countries_rand380.png")
plot(accum_country_380, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col="grey", main="Accumulation of GPSCs sampling 12 countries", xlab="No. countries", ylab="No. GPSCs",ylim=c(0,400))
dev.off()


#Select 12 random 380 samples from South Africa, storing GPSC counts per sample
#subset data to South Africa
ZA <- subset(gpsc_accum,Country=="South Africa", select = c(Taxon,GPSC))
#random shuffle of rows
ZA_set <- ZA[sample(nrow(ZA)),]
#split into samples of n=380
ZA_sets <- split(ZA_set, rep(1:ceiling(nrow(ZA)/380), each=380, length.out=nrow(ZA)))

#Create matrix for storing ZA results
random_sets <- data.frame(sort(unique(gpsc_accum$GPSC)))
row.names(random_sets) <- random_sets$sort.unique.gpsc_accum.GPSC..
random_sets[1] <- NULL

#Create ZA community dataset input for specaccum()
#merge ZA samples on GPSC rows to store GPSC counts per ZA sample set
i <- 1
for (sets in ZA_sets){
 format_380 <- as.data.frame(table(sets$GPSC))
 row.names(format_380) <- format_380$Var1
 format_380[1] <- NULL
 colnames(format_380) <- paste("set",i, "freq", sep=".")
 i <- i+1
 random_sets <- merge(random_sets,format_380,by="row.names",all=TRUE)
 rownames(random_sets) <- random_sets$Row.names
 random_sets$Row.names <- NULL
}
#replace NA with 0 and transpose
random_sets[is.na(random_sets)] <- 0
no_ZAGPSCs <- table(rowSums(random_sets)>0)
random_sets <- t(random_sets)

#run 100 random permutations of the ZA community dataset, where sites (random subsamples of ZA) are added in random order to determine GPSC accumulation curve
accum_ZA_380 <- specaccum(random_sets, method="random", permutations = 100)
png(file = "GPSC_accum_12samplesZA_rand380.png")
plot(accum_ZA_380, ci.type="poly", col="black", lwd=2, ci.lty=0, ci.col="grey", main="Accumulation of GPSCs using 12 random samples of South Africa", xlab="No. independant sample sets", ylab="No. GPSCs",ylim=c(0,400))
dev.off()


GPSC_test <- matrix(c(no_GPSCs[2],no_GPSCs[1],no_ZAGPSCs[2],no_ZAGPSCs[1]),nrow = 2)
#test if the number of GPSCs sampled in multiple countries is greater than a single location       
GPSC_sampled <- fisher.test(GPSC_test)$p    
