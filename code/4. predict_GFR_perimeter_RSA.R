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

# Script for making predictions from functional response models for burn perimeter RSA,
# as in Figure 4.

# Code written by Eric Palm

############################################################################################

require(glmmTMB)
require(dplyr)

# Load data for burn perimeter models
bp_2_week_summer <- readRDS("data/burn_perimeter_2_week.rds") %>% 
  filter(season == "summer")
bp_2_week_winter <- readRDS("data/burn_perimeter_2_week.rds") %>% 
  filter(season == "winter")

bp_24_hour_summer <- readRDS("data/burn_perimeter_24_hour.rds") %>% 
  filter(season == "summer")
bp_24_hour_winter <- readRDS("data/burn_perimeter_24_hour.rds") %>% 
  filter(season == "winter")

# Load functional response models for predictions
m_bp_2_week_summer <- readRDS("models/model_burn_perimeter_2_week_summer.rds")
m_bp_2_week_winter <- readRDS("models/model_burn_perimeter_2_week_winter.rds")
m_bp_24_hour_summer <- readRDS("models/model_burn_perimeter_24_hour_summer.rds")
m_bp_24_hour_winter <- readRDS("models/model_burn_perimeter_24_hour_winter.rds")

# Quick function to make sequence between min and max values of a covariate
make_seq <- function (x) {seq(min(x), max(x), len=100)}

# Create sample data frame for predictions across ranges of burn availability
# for each season-scale combination.
avail_s_burn_2_week <- 
  tibble(burn_avail_id_unscaled = make_seq(bp_2_week_summer$burn_avail_id),
         burn_avail_id = as.numeric(scale(burn_avail_id_unscaled, center = mean(bp_2_week_summer$burn_avail_id),
                               scale = sd(bp_2_week_summer$burn_avail_id))[,1]), 
         Season="Summer", Scale = "Two-week scale",
         TPI = 0, TRI = 0, tree = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0,
         lc_grass = 0, lc_fen = 0, lc_water = 0, lc_burn = 1, log_step_length = 0,
         herd = NA, id = NA, stratum = NA)

avail_w_burn_2_week <- 
  tibble(burn_avail_id_unscaled = make_seq(bp_2_week_winter$burn_avail_id),
         burn_avail_id = as.numeric(scale(burn_avail_id_unscaled, center = mean(bp_2_week_winter$burn_avail_id),
                                          scale = sd(bp_2_week_winter$burn_avail_id))[,1]), 
         Season="Winter", Scale = "Two-week scale",
         TPI = 0, TRI = 0, tree = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0,
         lc_grass = 0, lc_fen = 0, lc_water = 0, lc_burn = 1, log_step_length = 0,
         herd = NA, id = NA, stratum = NA)

avail_s_burn_24_hour <- 
  tibble(burn_avail_id_unscaled = make_seq(bp_24_hour_summer$burn_avail_id),
         burn_avail_id = as.numeric(scale(burn_avail_id_unscaled, center = mean(bp_24_hour_summer$burn_avail_id),
                                          scale = sd(bp_24_hour_summer$burn_avail_id))[,1]), 
         Season="Summer", Scale = "24-hour scale",
         TPI = 0, TRI = 0, tree = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0,
         lc_grass = 0, lc_fen = 0, lc_water = 0, lc_burn = 1, log_step_length = 0,
         herd = NA, id = NA, stratum = NA)

avail_w_burn_24_hour <- 
  tibble(burn_avail_id_unscaled = make_seq(bp_24_hour_winter$burn_avail_id),
         burn_avail_id = as.numeric(scale(burn_avail_id_unscaled, center = mean(bp_24_hour_winter$burn_avail_id),
                                          scale = sd(bp_24_hour_winter$burn_avail_id))[,1]), 
         Season="Winter", Scale = "24-hour scale",
         TPI = 0, TRI = 0, tree = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0,
         lc_grass = 0, lc_fen = 0, lc_water = 0, lc_burn = 1, log_step_length = 0,
         herd = NA, id = NA, stratum = NA)

# Function to predict from model using fixed effects only without intercept
pred_CI <- function(model, newdata=NULL, alpha=0.05) {
  X <- model.matrix(formula(model, fixed.only=TRUE)[-2], newdata)[,-1]
  pred0 <- as.vector(X %*% fixef(model)$cond[-1])
  V <- vcov(model)$cond[-1,-1]     
  pred_se <- sqrt(rowSums((X %*% V) * X)) 
  crit <- -qnorm(alpha/2)
  as_tibble(exp(cbind(conf_low=pred0-crit*pred_se,
                      conf_high=pred0+crit*pred_se, predict=pred0)))
}

# Make predictions for each data frame created above
preds_summer_2_week <- pred_CI(m_bp_2_week_summer, avail_s_burn_2_week)
preds_winter_2_week <- pred_CI(m_bp_2_week_winter, avail_w_burn_2_week)
preds_summer_24_hour <- pred_CI(m_bp_24_hour_summer, avail_s_burn_24_hour)
preds_winter_24_hour <- pred_CI(m_bp_24_hour_winter, avail_w_burn_24_hour)

# Join predictions with data frames for plotting and normalize between 0 and 1
all_preds <- bind_rows(
  bind_cols(avail_s_burn_2_week, preds_summer_2_week),
  bind_cols(avail_w_burn_2_week, preds_winter_2_week),
  bind_cols(avail_s_burn_24_hour, preds_summer_24_hour),
  bind_cols(avail_w_burn_24_hour, preds_winter_24_hour)) %>% 
  mutate(across(c(conf_low:predict), ~. /max(predict)))
