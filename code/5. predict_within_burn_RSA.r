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

# Script for making predictions from within-burn RSA, as in Figure 6.
# Code written by Eric Palm

############################################################################################

require(glmmTMB)
require(dplyr)

# Load in within-burn data 
wb_summer <- readRDS("data/within_burn.rds") %>% 
  filter(season == "summer")
wb_winter <- readRDS("data/within_burn.rds") %>% 
  filter(season == "winter")

# Load models for predictions
m_wb_summer <- readRDS("models/model_within_burn_summer.rds")
m_wb_winter <- readRDS("models/model_within_burn_winter.rds")

# Function to predict from glmmTMB model and drop intercept.
pred_CI <- function(model, newdata=NULL, alpha=0.05) {
  X <- model.matrix(formula(model, fixed.only=TRUE)[-2], newdata)[,-1]
  pred0 <- as.vector(X %*% fixef(model)$cond[-1])
  V <- vcov(model)$cond[-1,-1]     
  pred_se <- sqrt(rowSums((X %*% V) * X)) 
  crit <- -qnorm(alpha/2)
  as_tibble(exp(cbind(conf_low=pred0-crit*pred_se,
                      conf_high=pred0+crit*pred_se, predict=pred0)))
}

# Quick function to make sequence between min and max values of a covariate
make_seq <- function (x) {seq(min(x), max(x), len=100)}

#############################################################################
# Create sample data frame for predictions across ranges of covariates
# for each season-covariate combination.
#############################################################################

# Lichen

# Winter
lichen_predict_winter <- 
  tibble(lichen_unscaled = make_seq(wb_winter$lichen),
         lichen = as.numeric(scale(lichen_unscaled, center = mean(wb_winter$lichen),
                          scale = sd(wb_winter$lichen))[,1]), 
         TPI = 0, TRI = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0, lc_evergreen = 1, lc_grass = 0,
         lc_fen = 0, severity = 0, d_edge = 0, tsf = 0, id = NA, severity_01 = 1, stratum = NA)

# Summer
lichen_predict_summer <- 
  tibble(lichen_unscaled = make_seq(wb_summer$lichen),
         lichen = as.numeric(scale(lichen_unscaled, center = mean(wb_summer$lichen),
                          scale = sd(wb_summer$lichen))[,1]), 
         TPI = 0, TRI = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0, lc_evergreen = 1, lc_grass = 0, 
         lc_fen = 0, severity = 0, d_edge = 0, tsf = 0, id = NA, severity_01 = 1, stratum = NA)

# Make predictions
lichen_winter <- pred_CI(m_wb_winter, lichen_predict_winter)
lichen_summer <- pred_CI(m_wb_summer, lichen_predict_summer)

# Add season
fixed_winter_lichen <- bind_cols(lichen_predict_winter, as_tibble(lichen_winter)) %>% 
  mutate(x_range = lichen_unscaled,
         season = "Winter")

fixed_summer_lichen <- bind_cols(lichen_predict_summer, as_tibble(lichen_summer)) %>% 
  mutate(x_range = lichen_unscaled,
         season = "Summer")

# Limit prediction range for lichen to maximum 30% cover because values above this
# are very rare.

fixed_lichen <- 
  bind_rows(fixed_winter_lichen, fixed_summer_lichen) %>% 
  filter(x_range <= 30.1)
###############################################################

# Burn severity

# Winter
severity_predict_winter <- 
  tibble(severity_unscaled = make_seq(wb_winter$severity),
         severity = as.numeric(scale(severity_unscaled, center = mean(wb_winter$severity),
                                   scale = sd(wb_winter$severity))[,1]), 
         TPI = 0, TRI = 0, lichen = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0, lc_evergreen = 1,
         lc_grass = 0, lc_fen = 0, d_edge = 0, tsf = 0, id = NA, severity_01 = 1, stratum = NA)

# Summer
severity_predict_summer <- 
  tibble(severity_unscaled = make_seq(wb_summer$severity),
         severity = as.numeric(scale(severity_unscaled, center = mean(wb_summer$severity),
                                   scale = sd(wb_summer$severity))[,1]), 
         TPI = 0, TRI = 0, lichen = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0, lc_evergreen = 1,
         lc_grass = 0, lc_fen = 0, d_edge = 0, tsf = 0, id = NA, severity_01 = 1, stratum = NA)

# Make predictions
severity_winter <- pred_CI(m_wb_winter, severity_predict_winter)
severity_summer <- pred_CI(m_wb_summer, severity_predict_summer)

# Add season
fixed_winter_severity <- bind_cols(severity_predict_winter, as_tibble(severity_winter)) %>% 
  mutate(x_range = severity_unscaled,
         season = "Winter")

fixed_summer_severity <- bind_cols(severity_predict_summer, as_tibble(severity_summer)) %>% 
  mutate(x_range = severity_unscaled,
         season = "Summer")

# Limit prediction range for severity
fixed_severity <- 
  bind_rows(fixed_winter_severity, fixed_summer_severity) %>% 
  filter(x_range >= -100 & x_range <= 1000)
##############################################################

# Distance within burn perimeter

# Winter
d_edge_predict_winter <- 
  tibble(d_edge_unscaled = make_seq(wb_winter$d_edge),
         d_edge = as.numeric(scale(d_edge_unscaled, center = mean(wb_winter$d_edge),
                                   scale = sd(wb_winter$d_edge))[,1]),
         TPI = 0, TRI = 0, lichen = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0, lc_evergreen = 1, 
         lc_grass = 0, lc_fen = 0, severity = 0, tsf = 0, severity_01 = 1, id = NA, stratum = NA)

# Summer
d_edge_predict_summer <- 
  tibble(d_edge_unscaled = make_seq(wb_summer$d_edge),
         d_edge = as.numeric(scale(d_edge_unscaled, center = mean(wb_summer$d_edge),
                                   scale = sd(wb_summer$d_edge))[,1]),  
         TPI = 0, TRI = 0, lichen = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0, lc_evergreen = 1, 
         lc_grass = 0, lc_fen = 0, severity = 0, tsf = 0, severity_01 = 1, id = NA, stratum = NA)

# Make predictions
d_edge_winter <- pred_CI(m_wb_winter, d_edge_predict_winter)
d_edge_summer <- pred_CI(m_wb_summer, d_edge_predict_summer)

# Add season and convert units to kilometers
fixed_winter_d_edge <- bind_cols(d_edge_predict_winter, as_tibble(d_edge_winter)) %>% 
  mutate(x_range = d_edge_unscaled/1000,
         season = "Winter")

fixed_summer_d_edge <- bind_cols(d_edge_predict_summer, as_tibble(d_edge_summer)) %>% 
  mutate(x_range = d_edge_unscaled/1000,
         season = "Summer")

# Limit predictions range for distance within burn perimeter
fixed_d_edge <- 
  bind_rows(fixed_winter_d_edge, fixed_summer_d_edge) %>% 
  filter(x_range <= 7.6)

# Combine predictions from lichen, distance within burn perimeter, and severity.
# Normalize predictions between 0 and 1 for plotting.
top_panel_preds <-
  bind_rows(fixed_lichen, fixed_d_edge, fixed_severity) %>% 
  mutate(across(c(conf_low:predict), ~. /max(predict)))

#################################################################################

# Time since fire
tsf_predict_winter <- 
  tibble(tsf_unscaled = rep(make_seq(wb_winter$tsf), times = 3),
         tsf = as.numeric(scale(tsf_unscaled, center = mean(wb_winter$tsf),
                                scale = sd(wb_winter$tsf))[,1]), 
         TPI = 0, TRI = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0, lc_evergreen = 1,
         lc_grass = 0, lc_fen = 0, lichen = 0, d_edge = 0, severity_01 = 1,
         severity_unscaled = rep(c(0, 270, 900), each=100),
         severity = as.numeric(scale(severity_unscaled, 
                                     center = mean(wb_winter$severity),
                                     scale = sd(wb_winter$severity))[,1]), 
         id = NA, stratum = NA)

tsf_predict_summer <- 
  tibble(tsf_unscaled = rep(make_seq(wb_summer$tsf), times = 3),
         tsf = as.numeric(scale(tsf_unscaled, center = mean(wb_summer$tsf),
                                scale = sd(wb_summer$tsf))[,1]), 
         TPI = 0, TRI = 0, lc_shrubs = 0, lc_sparse = 0, lc_other = 0, lc_evergreen = 1,
         lc_grass = 0, lc_fen = 0, lichen = 0, d_edge = 0, severity_01 = 1,
         severity_unscaled = rep(c(0, 270, 900), each=100),
         severity = as.numeric(scale(severity_unscaled, 
                                     center = mean(wb_summer$severity),
                                     scale = sd(wb_summer$severity))[,1]), 
         id = NA, stratum = NA)

# Make predictions
tsf_winter <- pred_CI(m_wb_winter, tsf_predict_winter)
tsf_summer <- pred_CI(m_wb_summer, tsf_predict_summer)

# Add season
fixed_winter_tsf <- bind_cols(tsf_predict_winter, as_tibble(tsf_winter)) %>% 
  mutate(tsf_range = tsf_unscaled,
         season = "Winter")

fixed_summer_tsf <- bind_cols(tsf_predict_summer, as_tibble(tsf_summer)) %>% 
  mutate(tsf_range = tsf_unscaled,
         season = "Summer")

# Join together and add severity classes for plotting.
# Normalize predictions between 0 and 1 for plotting.
bottom_panel_preds <- 
  bind_rows(fixed_winter_tsf, fixed_summer_tsf) %>% 
  mutate(severity_class = if_else(severity < -1, "Unburned residuals",
                                if_else(severity > 1, "High severity", "Low severity")),
         severity_class = factor(severity_class, 
                                 levels = c("Unburned residuals", "Low severity", "High severity"))) %>% 
  mutate(across(c(conf_low:predict), ~. /max(predict)))
