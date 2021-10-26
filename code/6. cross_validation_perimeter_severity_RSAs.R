
############################################################################################
#                                                                                           
# Journal: Ecological Applications
#
# Article title: "Increasing fire frequency and severity will increase habitat loss 
#                 for a boreal forest indicator species"
#
# Authors: Eric C. Palm, Michael J. Suitor, Kyle C. Joly, Jim D. Herriges, Allicia P. Kelly, 
#          Dave Hervieux, Kelsey L.M. Russell, Torsten W. Bentzen, Nicholas C. Larter, 
#          and Mark Hebblewhite
#
############################################################################################

# Script for running cross-validation for burn perimeter and burn severity RSAs
# Code written by Eric Palm

# Spearman rank correlations appear in Appendix S1 Table S8 and adjusted frequencies
# of used locations within SSF bins appear in Appendix S1 Figure S5.

############################################################################################

require(caret)
require(glmmTMB)
require(tidyverse)

# Change these depending on the desired analysis
m_season <- "summer" # or "winter
m_scale <- "24_hour" # or "2_week"
m_type <- "burn_severity" # or "burn_perimeter"

# Function to predict mean values from glmmTMB model, dropping intercept.
pred_fun <- function(model, newdata=NULL) {
  X <- model.matrix(formula(model, fixed.only=TRUE)[-2], newdata)[,-1]
  exp(as.vector(X %*% fixef(model)$cond[-1]))
}

# Load in SSF data frame.
df <- readRDS(paste0("data/", m_type, "_", m_scale, ".rds")) %>% 
  filter(season == m_season) 

# Load in top model.
model_CV <- readRDS(paste0("models/model_", m_type, "_", m_scale, "_", m_season, ".rds"))

# Get model formula
model_formula <- formula(model_CV)

# Create training data folds, each of which is 14 of 15 populations (herds)
# and scale all numeric (double) variables.
train_list <- list()
for (i in levels(as.factor(df$herd))) {
  
  train_list[[i]] <- filter(df, herd != i) %>% 
    mutate(across(where(is.double), ~as.numeric(scale(.)[,1])))
  
}

# Create test data folds, each of which is the 1 population withheld from
# the corresponding training fold.
# Apply the scaling from the training data to the corresponding test data.
test_list <- list()
for (i in levels(as.factor(df$herd))) {
  
  others <- df %>% 
    filter(herd == i) %>% 
    dplyr::select_if(negate(is.double))
  
  train_params <- df %>% 
    filter(herd != i) %>% 
    dplyr::select_if(is.double) %>% 
    caret::preProcess(., method = c("center", "scale"))
  
  test_list[[i]] <- df %>% 
    filter(herd == i) %>% 
    dplyr::select_if(is.double) %>% 
    predict(train_params, .) %>% 
    bind_cols(others, .)
  
}

# Create data frame from these lists where each fold is the withheld population.
folds <- tibble(fold=levels(as.factor(df$herd)),
                train_data= train_list,
                test_data=test_list)

# Fit the model to each training fold and save it in the data frame.
system.time(
  folds_models <- folds %>%
    mutate(model = purrr::map(train_data, ~glmmTMB(model_formula, data=., family=poisson,
                                            map=list(theta=factor(c(NA, 1:nlevels(model_CV$modelInfo$map$theta)))),
                                            start=list(theta=c(log(1e3), rep(0, nlevels(model_CV$modelInfo$map$theta))))))) 
)                 

# Use the model fit from the training data to predict over the test data
# and create a data frame with just these predictions.
system.time(
  folds_predict <-
    folds_models %>% 
    mutate(predict = purrr::map2(model, test_data, ~pred_fun(.x, newdata=.y))) %>% 
    dplyr::select(fold, test_data, predict) %>% 
    unnest(cols=c(test_data, predict))
)

# Rank each test fold prediction from 1-11 by stratum (1 used, 10 available locations).
# Calculate adjusted frequencies of used locations within each bin.
binned_preds <-
  folds_predict %>% 
  group_by(fold, stratum) %>% 
  mutate(bin = ntile(predict, 11)) %>% 
  group_by(fold, bin, case_) %>% 
  dplyr::summarize(count=n()) %>% 
  group_by(fold, case_) %>% 
  mutate(prop=count/sum(count)) %>% 
  group_by(fold, bin) %>% 
  dplyr::summarize(adjusted_freq = prop[case_==T]/prop[case_==F]) %>% 
  ungroup() 

# Calculate Spearman rank correlation for each population
binned_preds %>% 
  group_by(Population = fold) %>% 
  dplyr::summarize(spearman = cor.test(bin, adjusted_freq, method="spearman")$estimate)

# Save results 
saveRDS(binned_preds, file=paste0("results/", m_type, "_", m_scale, "_", m_season, "_binned_CV.rds"))


