paper2 <- read.csv("Paper2-supplementary.csv", header = TRUE, sep =",")

countries <- c("Malawi", "The Gambia", "China")

GPSC_NVT_prop <- matrix(data=NA,nrow=0,ncol=8)
for (x in countries){
  #Determine GPSC type per country based on VT proportion (<50% or >=50%) in first period observed
  GPSC_props <- matrix(data=NA,nrow=0,ncol=3)
  colnames(GPSC_props) <- c("GPSC","Proportion VT", "type")
  coun_p2 <- subset(paper2, Country==x, select=c(GPSC,Vaccine_Period,Vaccine_Status,In_Silico_Serotype,Country,Year))
  GPSCs <- unique(coun_p2$GPSC)
  for (GPSC in GPSCs){
    GPSC_data <- subset(coun_p2, GPSC==GPSC)
    #table of vaccine period by vaccine status, remove rows that sum to sero (non isolates from that period)
    GPSC_tab <- as.data.frame(unclass(table(GPSC_data$Vaccine_Period,GPSC_data$Vaccine_Status)))[ rowSums(as.data.frame(unclass(table(GPSC_data$Vaccine_Period,GPSC_data$Vaccine_Status))))!=0, ]
    #Calcuate the proportion VT in the earliest period found in last row of table
    GPSC_prop <- tail(GPSC_tab$PCV/rowSums(GPSC_tab), n=1)
    GPSC_props <- rbind(GPSC_props, c(GPSC,GPSC_prop,ifelse(GPSC_prop<0.5,"NVT","VT")))
  }
  GPSC_props <- as.data.frame(GPSC_props)
  `NVT-GPSCs` <- as.character(subset(GPSC_props, type=="NVT", select=c(GPSC))$GPSC)
  `VT-GPSCs` <- as.character(subset(GPSC_props, type=="VT", select=c(GPSC))$GPSC)
  GPSC_type <- list(`NVT-GPSCs`,`VT-GPSCs`)
  names(GPSC_type) <- c("NVT-GPSC","VT-GPSC")
  
  dat <- subset(paper2, Country==x & Vaccine_Period!="Post-PCV7" & Vaccine_Status_extended=="NVT", select = c(GPSC, Vaccine_Period))
  GPSCs <- as.data.frame.matrix(table(droplevels(dat)))
  colsum <- colSums(GPSCs)
  GPSCs$tot <- GPSCs$`Post-PCV13`+GPSCs$`Pre-PCV`
  GPSCs <- subset(GPSCs, tot>3, select = c(`Post-PCV13`,`Pre-PCV`))
  if(nrow(GPSCs) >0){
    for (i in 1:nrow(GPSCs)){
     row <- GPSCs[i,]
     by2 <- matrix(c(row[1,2],colsum[2]-row[1,2],row[1,1],colsum[1]-row[1,1]),ncol=2,byrow=TRUE)
     p <- fisher.test(by2, alternative)$p.value
     #what type of GPSC is being tested?
     if (rownames(row) %in% GPSC_type$`NVT-GPSC`){
       what_type <- "NVT-GPSC"
     } else {
       what_type <- "VT-GPSC"
     }
     GPSC_NVT_prop <- rbind(GPSC_NVT_prop,c(x, rownames(row), what_type, row[1,2],colsum[2],row[1,1],colsum[1],p))
    }
  }
}
#n=10 tests no multiple testing
colnames(GPSC_NVT_prop) <- c("Country","GPSC","GPSC_type","Pre","Pre total","PCV13","PCV13 total","p-value")
write.csv(GPSC_NVT_prop, file ="Tab3_GPSC_NVT_postprop.csv", row.names = FALSE)
