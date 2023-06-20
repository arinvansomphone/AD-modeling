# loading packages 
library(tidyverse)
library(readr)
library(caTools)
library(caret)
library(glmnet)
library(faux)

# importing and tidying the data
upennbiomarkers <- read_csv("R work/dementia modeling/UPENNPLASMA_29May2023.csv")
ADNItraining1 <- read_csv("R work/dementia modeling/AD_Challenge_Training_Data_Clinical_Updated_7.22.2014/ADNI_Training_Q3_APOE_CollectionADNI1Complete 1Yr 1.5T_July22.2014.csv")
ADNItraining2 <- read_csv("R work/dementia modeling/AD_Challenge_Training_Data_Clinical_Updated_7.22.2014/ADNI_Training_Q1_APOE_July22.2014.csv")
blennowplasmatau <- read_csv("R work/dementia modeling/BLENNOWPLASMATAU_29May2023.csv")
RADERlipidomics <- read_csv("R work/dementia modeling/ADNI_LIPIDOMICSRADER_29May2023.csv") 
DESIKANLABPHS <- read_csv("R work/dementia modeling/DESIKANLAB_03Jun2023.csv")
homocys <- read_csv("R work/dementia modeling/HCRES_03Jun2023.csv")

upennbiomarkers <- upennbiomarkers %>%
  filter(VISCODE == "bl") %>%
  select(RID, AB40, AB42)

blennowplasmatau <- blennowplasmatau %>%
  select(RID, PLASMATAU)

RADERlipidomics <- RADERlipidomics %>%
  select(RID, CHOL, HDL, TG, APOA1, APOE)

ADNItraining1 <- ADNItraining1 %>%
  select(RID, DX.bl, AGE, PTGENDER, PTEDUCAT, PTETHCAT, PTRACCAT,
         APOE4, MMSE, `APOE Genotype`) %>%
  rename(diagnosis = DX.bl)

ADNItraining2 <- ADNItraining2 %>%
  filter(VISCODE == "m24") %>%
  select(RID, MMSE) %>%
  rename(MMSE24 = MMSE)

DESIKANLABPHS <- DESIKANLABPHS %>%
  select(RID, PHS, CIR)

homocys <- homocys %>%
  filter(Phase == "ADNI1", VISCODE == "bl") %>%
  select(RID, HCAMPLAS, HCVOLUME) 

completedata <- left_join(ADNItraining1, upennbiomarkers, 
                          by = join_by(RID == RID))
completedata <- left_join(completedata, blennowplasmatau, 
                          by = join_by(RID == RID))
completedata <- left_join(completedata, RADERlipidomics,
                          by = join_by(RID == RID))
completedata <- left_join(completedata, DESIKANLABPHS,
                          by = join_by(RID ==RID))
completedata <- left_join(completedata, homocys,
                          by = join_by(RID ==RID))
completedata <- left_join(completedata, ADNItraining2,
                          by = join_by(RID == RID))

completedata <- completedata %>% 
  select(-RID) %>%
  rename(APOEgenotype = `APOE Genotype`)
names <- c(1, 3, 5, 6, 7, 9)
completedata[,names] <- lapply(completedata[,names], factor)
completedata$CHOL <- as.numeric(completedata$CHOL)
completedata$HDL <- as.numeric(completedata$HDL)
completedata$TG <- as.numeric(completedata$TG)
completedata$APOA1 <- as.numeric(completedata$APOA1)
completedata$APOE <- as.numeric(completedata$APOE)
completedata <- completedata %>%
  na.omit()

# logistic regression model 
## more tidying and making variables factors
completedata$diagnosis <- str_replace(completedata$diagnosis, "LMCI", "AD")
completedata$diagnosis <- as.factor(completedata$diagnosis)

## train test set creation 
set.seed(100)
spl = sample.split(completedata$diagnosis, SplitRatio = 0.7)
train = subset(completedata, spl == TRUE)
test = subset(completedata, spl == FALSE)

## model creation 
model_glm = glm(diagnosis ~ ., family = "binomial", data = train)
summary(model_glm)

## predictions on the training set
predictTrain = predict(model_glm, data = train, type = "response")
train_and_predictions <- train %>%
  mutate(predicttrain = predictTrain)

## confusion matrix on training data
table(train$diagnosis, predictTrain >= 0.5)
(165+89)/nrow(train) ### 87% accuracy 

## predictions on the test set
predictTest = predict(model_glm, newdata = test, type = "response")

## confusion matrix on the test set
table(test$diagnosis, predictTest >= 0.5)
(70+38)/nrow(test) ### 83% accuracy

# decision tree
## making train and test datasets
set.seed(245)
spl = sample.split(completedata$diagnosis, SplitRatio = 0.7)
train = subset(completedata, spl == TRUE)
test = subset(completedata, spl == FALSE)

## training the network
library(rpart)
library(rpart.plot)

fit <- rpart(diagnosis ~ ., data = train, method = 'class')
rpart.plot(fit, extra = 106)

predict_unseen <- predict(fit, test, type = 'class')

table_mat <- table(test$diagnosis, predict_unseen)
table_mat
sum(diag(table_mat)) / sum(table_mat) ### 75% accuracy

## tuning
accuracy_tune <- function(fit) {
  predict_unseen <- predict(fit, test, type = 'class')
  table_mat <- table(test$diagnosis, predict_unseen)
  accuracy_Test <- sum(diag(table_mat)) / sum(table_mat)
  accuracy_Test
}
control <- rpart.control(minsplit = 4, minbucket = round(5 / 3), maxdepth = 3,
                         cp = 0)
tune_fit <- rpart(diagnosis ~ ., data = train, method = 'class',
                  control = control)
accuracy_tune(tune_fit) ### 67% accuracy, worsens with tuning 

# multivariate multiple regression model
## simulating more data 
simulateddata <- sim_df(completedata, 1000)

## removing categorical variables and joining with actual data
simulateddataregression <- completedata %>%
  select(-c(diagnosis, PTRACCAT, PTGENDER, PTETHCAT, APOEgenotype, APOE4))
simulateddata <- simulateddata %>%
  select(-c(id))
totalsimulated <- rbind(simulateddataregression, simulateddata)

## select factors desired
completedataregression <- completedata %>%
  select(-c(diagnosis, PTRACCAT)) %>%
  filter(!APOEgenotype == 22)

## making train and test datasets
set.seed(245)
spl = sample.split(totalsimulated$MMSE24, SplitRatio = 0.7)
train = subset(totalsimulated, spl == TRUE)
test = subset(totalsimulated, spl == FALSE)

mlml <- lm(MMSE24 ~ ., data = train)
summary(mlml)$coefficient
sigma(mlml) / mean(train$MMSE24)

## RMSE = root mean square error
simulatedregression_rmse_train <- train %>%
  mutate(pred.MMSE24.train = predict(mlml))

plot(simulatedregression_rmse_train$MMSE24, 
     simulatedregression_rmse_train$pred.MMSE24.train)

mse <- simulatedregression_rmse_train %>%
  mutate(error = pred.MMSE24.train - MMSE24,
         sq.error = error^2) %>%
  summarise(mse = mean(sq.error))
rmse <- sqrt(mse) ### training rmse = 3.50
rmse

## check rmse on test data
pred_reg_test <- predict(mlml, newdata = test)
reg_rmse_test <- sqrt(mean((pred_reg_test - test$MMSE24)^2))
reg_rmse_test ### test rmse = 3.43

## LASSO model
y <- train$MMSE24
x <- data.matrix(train[, 1:15])
cv_model <- cv.glmnet(x, y, alpha = 1)

best_lambda <- cv_model$lambda.min
best_lambda

plot(cv_model)

best_model <- glmnet(x, y, alpha = 1, lambda = best_lambda)
coef(best_model)

simulatedlasso_rmse_train <- train %>%
  mutate(pred.MMSE24.train = predict(best_model, s = best_lambda, 
                                     newx = x))

## RMSE of LASSO model
plot(simulatedlasso_rmse_train$MMSE24, 
     simulatedlasso_rmse_train$pred.MMSE24.train)

lassomse <- simulatedlasso_rmse_train %>%
  mutate(error = pred.MMSE24.train - MMSE24,
         sq.error = error^2) %>%
  summarise(mse = mean(sq.error))
lassormse <- sqrt(lassomse) ### training rmse = 3.42

## check RMSE on test data 
z <- data.matrix(test[, 1:15])

pred_lasso_test <- predict(best_model, s = best_lambda, newx = z)
lasso_rmse_test <- sqrt(mean((pred_lasso_test - test$MMSE24)^2))
lasso_rmse_test ### test rmse = 3.23
