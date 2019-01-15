#input data T8-GPSC_summaries tab from Supplementary-T1-T21 in csv format
T8 <- read.csv("T8-GPSC_summaries.csv", header = TRUE, sep =",")
T8 <- subset(T8, select = c("GPSC", "Mean_No_resistant_classes"))
#input data T11-GPSC_summaries tab from Supplementary-T1-T21 in csv format
T11 <- read.csv("T11-GPSC_recombination.csv", header = TRUE, sep =",")
T11 <- subset(T11, select = c("GPSC", "rho.theta","r.m"))

#combine input
data <- merge(x=subset(T8, GPSC %in% unique(T11$GPSC)),y=T11, by = "GPSC")

#transform dependant
classes <- sqrt(data$Mean_No_resistant_classes)

#linear regression to predict continous mean number of resistant antibiotic classes
lmrho <- lm(classes ~ data$rho.theta)
lmrm <- lm(classes ~ data$r.m)

#is there heteroskedasticity (unequal variances of the residuals)
par(mfrow=c(2,2))
plot(lmrho2)
lmtest::bptest(lmrho)
car::ncvTest(lmrho)

#model
summary(lmrho)
summary(lmrm)
