#input isolate data from Supplementary data
paper2 <- read.csv("Paper2-supplementary_v4.csv", header = TRUE, sep =",")

#Pen
dat <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(dat$Vaccine_Period)[c(2,1)])
dat$Predicted_PEN_MIC_CLSI <- factor(dat$Predicted_PEN_MIC_CLSI,levels(dat$Predicted_PEN_MIC_CLSI)[c(2,1)])
glm_penR <- glm(Predicted_PEN_MIC_CLSI ~ Country + Vaccine_Period, data = dat, family = binomial)

datNVT <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7" & Vaccine_Status=="NVT"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(datNVT$Vaccine_Period)[c(2,1)])
datNVT$Predicted_PEN_MIC_CLSI <- factor(datNVT$Predicted_PEN_MIC_CLSI,levels(datNVT$Predicted_PEN_MIC_CLSI)[c(2,1)])
glm_penR_NVT <- glm(Predicted_PEN_MIC_CLSI ~ Vaccine_Period + Country, data = datNVT, family = binomial)

#Ery
dat <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(dat$Vaccine_Period)[c(2,1)])
dat$Predicted_Erythromycin_susceptibility <- factor(dat$Predicted_Erythromycin_susceptibility,levels(dat$Predicted_Erythromycin_susceptibility)[c(2,1)])
glm_EryR <- glm(Predicted_Erythromycin_susceptibility ~ Country + Vaccine_Period, data = dat, family = binomial)

datNVT <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7" & Vaccine_Status=="NVT"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(datNVT$Vaccine_Period)[c(2,1)])
datNVT$Predicted_Erythromycin_susceptibility <- factor(datNVT$Predicted_Erythromycin_susceptibility,levels(datNVT$Predicted_Erythromycin_susceptibility)[c(2,1)])
glm_EryR_NVT <- glm(Predicted_Erythromycin_susceptibility ~ Country + Vaccine_Period, data = datNVT, family = binomial)


#Tet
dat <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(dat$Vaccine_Period)[c(2,1)])
dat$Predicted_Tetracycline_susceptibility <- factor(dat$Predicted_Tetracycline_susceptibility,levels(dat$Predicted_Tetracycline_susceptibility)[c(2,1)])
glm_tetR <- glm(Predicted_Tetracycline_susceptibility ~ Country + Vaccine_Period, data = dat, family = binomial)


datNVT <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7" & Vaccine_Status=="NVT"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(datNVT$Vaccine_Period)[c(2,1)])
datNVT$Predicted_Tetracycline_susceptibility <- factor(datNVT$Predicted_Tetracycline_susceptibility,levels(datNVT$Predicted_Tetracycline_susceptibility)[c(2,1)])
glm_tetR_NVT <- glm(Predicted_Tetracycline_susceptibility ~ Vaccine_Period + Country, data = datNVT, family = binomial)

#SXT
dat <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7" ))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(dat$Vaccine_Period)[c(2,1)])
dat$Predicted_Cotrimoxazole_susceptibility <- factor(dat$Predicted_Cotrimoxazole_susceptibility,levels(dat$Predicted_Cotrimoxazole_susceptibility)[c(2,1)])
glm_SXT_NS <- glm(Predicted_Cotrimoxazole_susceptibility ~ Country + Vaccine_Period, data = dat, family = binomial)

datNVT <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7" & Vaccine_Status=="NVT"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(datNVT$Vaccine_Period)[c(2,1)])
datNVT$Predicted_Cotrimoxazole_susceptibility <- factor(datNVT$Predicted_Cotrimoxazole_susceptibility,levels(datNVT$Predicted_Cotrimoxazole_susceptibility)[c(2,1)])
glm_SXT_NS_NVT <- glm(Predicted_Cotrimoxazole_susceptibility ~ Vaccine_Period + Country, data = datNVT, family = binomial)

#Chl
dat <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(dat$Vaccine_Period)[c(2,1)])
dat$Predicted_Chloramphenicol_susceptibility <- factor(dat$Predicted_Chloramphenicol_susceptibility,levels(dat$Predicted_Chloramphenicol_susceptibility)[c(2,1)])
glm_ChlR <- glm(Predicted_Chloramphenicol_susceptibility ~ Country + Vaccine_Period, data = dat, family = binomial)

datNVT <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7" & Vaccine_Status=="NVT"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(datNVT$Vaccine_Period)[c(2,1)])
datNVT$Predicted_Chloramphenicol_susceptibility <- factor(datNVT$Predicted_Chloramphenicol_susceptibility,levels(datNVT$Predicted_Chloramphenicol_susceptibility)[c(2,1)])
glm_ChlR_NVT <- glm(Predicted_Chloramphenicol_susceptibility ~ Vaccine_Period + Country, data = datNVT, family = binomial)


#Number of classes
dat <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7"))
dat$Vaccine_Period <- factor(dat$Vaccine_Period,levels(dat$Vaccine_Period)[c(2,1)])
glm_MDR <- glm(MDR ~ Country + Vaccine_Period, data = dat, family = binomial)

datNVT <- droplevels(subset(paper2, Vaccine_Period!="Post-PCV7" & Vaccine_Status=="NVT"))
datNVT$Vaccine_Period <- factor(datNVT$Vaccine_Period,levels(datNVT$Vaccine_Period)[c(2,1)])
glm_MDR_NVT <- glm(MDR ~ Vaccine_Period + Country, data = datNVT, family = binomial)

results <-rbind(
c("penR_post",summary(glm_penR)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("ChlR_post",summary(glm_ChlR)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("EryR_post",summary(glm_EryR)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("SXTR_post",summary(glm_SXT_NS)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("tetR_post",summary(glm_tetR)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("MDR_post",summary(glm_MDR)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("penR_NVTpost",summary(glm_penR_NVT)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("ChlR_NVTpost",summary(glm_ChlR_NVT)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("EryR_NVTpost",summary(glm_EryR_NVT)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("SXTR_NVTpost",summary(glm_SXT_NS_NVT)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("tetR_NVTpost",summary(glm_tetR_NVT)$coefficients['Vaccine_PeriodPost-PCV13',]),
c("MDR_NVTpost",summary(glm_MDR_NVT)$coefficients['Vaccine_PeriodPost-PCV13',])
)

colnames(results)[1] <- "Type"
write.csv(results, file ="Tab4_post_NVT_resistance.csv", row.names = FALSE)
