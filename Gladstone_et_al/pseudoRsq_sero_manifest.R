#install.packages("rcompanion")
#install.packages("glmnet")
#install.packages("fmsb")

library("rcompanion")
library("glmnet")
library("fmsb")

#input T4 tab from Supplementary T1-T21 in csv format
T4 <- read.csv("T4-Heritability_dataset.csv", row.names=1, header = TRUE, sep =",", stringsAsFactors = FALSE)

sero_data <- subset(T4, select=c(Clinical.manifest,Serotype))
sero_data$Clinical.manifest[sero_data$Clinical.manifest == "Disease"] <- 1
sero_data$Clinical.manifest[sero_data$Clinical.manifest == "Carriage"] <- 0
sero_data$Clinical.manifest <- as.factor(sero_data$Clinical.manifest)
sero <- model.matrix(~Serotype, data=sero_data)[,-1]

#negelkerke pseudo r2 (logistic model) without selector prediction

model = glm(sero_data$Clinical.manifest ~ sero, family = 'binomial')

NagelkerkeR2(model)

#lasso regression with leave one out cross validation
lasso_X <- sero
lasso_y <- sero_data$Clinical.manifest
lasso.sero <- glmnet(lasso_X, lasso_y, alpha = 1,
                     family = 'binomial')
plot(lasso.sero,xvar="lambda",label=T)
cv.lasso.sero <- cv.glmnet(lasso_X, lasso_y, family = 'binomial',
                           alpha = 1, nfolds = length(lasso_y))
plot(cv.lasso.sero)
coef(cv.lasso.sero, s="lambda.1se")

# selected predictors in regression
selected <- lasso_X[,which(coef(cv.lasso.sero, s="lambda.1se")[-1] != 0)]

#negelkerke pseudo r2 (logistic model) with selector prediction
glm.sero <- glm(lasso_y ~ selected, family = binomial())

NagelkerkeR2(glm.sero)