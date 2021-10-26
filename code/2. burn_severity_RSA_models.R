
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

# Script for running top models in burn severity RSA
# Code written by Eric Palm

############################################################################################

require(glmmTMB)
require(dplyr)

# Load data
bs_24_hour <- readRDS("data/burn_severity_24_hour.rds") 

# Create seasonal datasets with centered and scaled continuous covariates
bs_24_hour_summer <-
  bs_24_hour %>% 
  filter(season == "summer") %>% 
  mutate(across(c(TRI:tree, log_step_length), ~as.numeric(scale(.)[,1])))

bs_24_hour_winter <-
  bs_24_hour %>% 
  filter(season == "winter") %>% 
  mutate(across(c(TRI:tree, log_step_length), ~as.numeric(scale(.)[,1])))

#######################################################################################
# Top burn severity models.

# Fixed-effects coefficients for burn severity models appear in Figure 5 and in
# Appendix S1 Table S4B.

# Random coefficients for burn severity categories appear in Appendix S1 Table S7.

# Note that 'map' and 'start' arguments tell glmmTMB to fix the variance of the random
# intercept at log(1e3) and freely estimate all nine random slopes in the models.
#######################################################################################

# 24 HOUR MODELS

m_bs_24_hour_summer <-
  glmmTMB(case_ ~ tree  + TRI + TPI + 
            lc_shrubs + lc_grass + lc_sparse + lc_regen +
            lc_residual + lc_old_burn + lc_low_severity +
            lc_high_severity +
            lc_water + lc_other + lc_fen + log_step_length + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_regen | herd/id) + (0 + lc_residual | herd/id) + 
            (0 + lc_old_burn | herd/id) + (0 + lc_low_severity | herd/id) + 
            (0 + lc_high_severity | herd/id),
          data=bs_24_hour_summer, family=poisson,
          map=list(theta=factor(c(NA, 1:17))),
          start=list(theta=c(log(1e3), rep(0, 17))))


m_bs_24_hour_winter <-
  glmmTMB(case_ ~ tree + I(tree^2) + TRI + I(TRI^2) + TPI + 
            lc_shrubs + lc_grass + lc_sparse + lc_regen +
            lc_residual + lc_old_burn + lc_low_severity +
            lc_high_severity +
            lc_water + lc_other + lc_fen + log_step_length + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_regen | herd/id) + (0 + lc_residual | herd/id) + 
            (0 + lc_old_burn | herd/id) + (0 + lc_low_severity | herd/id) + 
            (0 + lc_high_severity | herd/id),
          data=bs_24_hour_winter, family=poisson,
          map=list(theta=factor(c(NA, 1:17))),
          start=list(theta=c(log(1e3), rep(0, 17))))

saveRDS(m_bs_24_hour_summer, "models/model_burn_severity_24_hour_summer.rds")
saveRDS(m_bs_24_hour_winter, "models/model_burn_severity_24_hour_winter.rds")


