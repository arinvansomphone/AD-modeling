# loading packages
library(glmnet)
library(randomForest)
library(ranger)
library(tidymodels)
library(readr)

# importing and tidying the data
upennbiomarkers <- read_csv("R work/dementia modeling/UPENNPLASMA_29May2023.csv")
ADNItraining1 <- read_csv("R work/dementia modeling/AD_Challenge_Training_Data_Clinical_Updated_7.22.2014/ADNI_Training_Q3_APOE_CollectionADNI1Complete 1Yr 1.5T_July22.2014.csv")
ADNItraining2 <- read_csv("R work/dementia modeling/AD_Challenge_Training_Data_Clinical_Updated_7.22.2014/ADNI_Training_Q1_APOE_July22.2014.csv")
blennowplasmatau <- read_csv("R work/dementia modeling/BLENNOWPLASMATAU_29May2023.csv")
RADERlipidomics <- read_csv("R work/dementia modeling/ADNI_LIPIDOMICSRADER_29May2023.csv") 
DESIKANLABPHS <- read_csv("R work/dementia modeling/DESIKANLAB_03Jun2023.csv")
homocys <- read_csv("R work/dementia modeling/HCRES_03Jun2023.csv")
baselinesymptoms <- read_csv("R work/dementia modeling/BLSCHECK_08Jul2023.csv")
neuroexam <- read_csv("R work/dementia modeling/NEUROEXM_08JuL2023.csv")

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

baselinesymptoms <- baselinesymptoms %>%
  filter(Phase == "ADNI1") %>%
  select(RID, BCNAUSEA, BCVOMIT, BCDIARRH, BCCONSTP, BCABDOMN, BCSWEATN,
         BCDIZZY, BCENERGY, BCDROWSY, BCVISION, BCHDACHE, BCDRYMTH, BCBREATH,
         BCCOUGH, BCPALPIT, BCCHEST, BCURNDIS, BCANKLE, BCMUSCLE, BCRASH, 
         BCINSOMN, BCDPMOOD, BCCRYING, BCELMOOD, BCWANDER, BCFALL, BCOTHER)

neuroexam <- neuroexam %>%
  filter(Phase == "ADNI1") %>%
  select(RID, NXVISUAL, NXAUDITO, NXTREMOR, NXCONSCI, NXNERVE, NXMOTOR, 
         NXFINGER, NXHEEL, NXSENSOR, NXTENDON, NXPLANTA, NXGAIT, NXOTHER,
         NXABNORM)

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
completedata <- left_join(completedata, baselinesymptoms,
                          by = join_by(RID ==RID))
completedata <- left_join(completedata, neuroexam,
                          by = join_by(RID ==RID))
completedata <- left_join(completedata, ADNItraining2,
                          by = join_by(RID == RID))

completedata <- completedata %>% 
  select(-c(RID, diagnosis, PTRACCAT, PTETHCAT)) %>% 
  rename(APOEgenotype = `APOE Genotype`) %>%
  # remove PTRACCAT AND PTETHCAT due to data issues 
  filter(!APOEgenotype == 22) %>%
  na.omit()
  # remove due to data issues
names <- c(1, 3, 5, 6, 7, 9)
completedata[,names] <- lapply(completedata[,names], factor)
completedata$CHOL <- as.numeric(completedata$CHOL)
completedata$HDL <- as.numeric(completedata$HDL)
completedata$TG <- as.numeric(completedata$TG)
completedata$APOA1 <- as.numeric(completedata$APOA1)
completedata$APOE <- as.numeric(completedata$APOE)

# splitting training and testing data
set.seed(4595)
data_split <- initial_split(completedata, strata = "MMSE24", prop = 0.75)

train <- training(data_split)
test <- testing(data_split)

# random forest model 
rf_defaults <- rand_forest(mode = "regression")

preds <- colnames(completedata)[1:59]

rf_xy_fit <- 
  rf_defaults %>%
  set_engine("ranger") %>%
  fit_xy(
    x = train[, preds],
    y = log10(train$MMSE24)
  )

test_results <-
  test %>%
  select(MMSE24) %>%
  mutate(MMSE24 = log10(MMSE24)) %>%
  bind_cols(
    predict(rf_xy_fit, new_data = test[, preds])
  )
test_results %>% slice(1:5)

## summarize performance 
test_results %>% metrics(truth = MMSE24, estimate = .pred)

rand_forest(mode = "regression", mtry = .preds(), trees = 1000) %>%
  set_engine("ranger") %>%
  fit(
    log10(MMSE24) ~ .,
    data = train
  )

# Regularization
norm_recipe <-
  recipe(
    MMSE24 ~ .,
    data = train
  ) %>%
  step_dummy(all_nominal()) %>%
  step_center(all_predictors()) %>%
  step_scale(all_predictors()) %>%
  step_log(MMSE24, base = 10) %>%
  prep(training = train, retain = TRUE)

## fit model using processed version of data
glmn_fit <- 
  linear_reg(penalty = 0.001, mixture = 0.5) %>%
  set_engine("glmnet") %>%
  fit(MMSE24 ~., data = bake(norm_recipe, new_data = NULL))
glmn_fit

# finding predictions for glmn
test_normalized <- bake(norm_recipe, new_data = test, all_predictors())

test_results <-
  test_results %>%
  rename(`random forest` = .pred) %>%
  bind_cols(
    predict(glmn_fit, new_data = test_normalized) %>%
      rename(glmnet = .pred)
  )
test_results

test_results %>% metrics(truth = MMSE24, estimate = glmnet)

test_results %>%
  gather(model, prediction, -MMSE24) %>%
  ggplot(aes(x = prediction, y = MMSE24)) +
  geom_abline(col = "green", lty = 2) +
  geom_point(alpha = 0.8) +
  facet_wrap(~model) +
  coord_fixed()

# https://www.tidymodels.org/start/resampling/