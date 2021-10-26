
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

# Script for running top models in burn perimeter RSA
# Code written by Eric Palm

############################################################################################

require(glmmTMB)
require(dplyr)

# Load data
bp_2_week <- readRDS("data/burn_perimeter_2_week.rds") 
bp_24_hour <- readRDS("data/burn_perimeter_24_hour.rds") 

# Create seasonal datasets with centered and scaled continuous covariates
bp_2_week_summer <-
  bp_2_week %>% 
  filter(season == "summer") %>% 
  mutate(across(where(is.double), ~as.numeric(scale(.)[,1])))

bp_2_week_winter <-
  bp_2_week %>% 
  filter(season == "winter") %>% 
  mutate(across(where(is.double), ~as.numeric(scale(.)[,1])))

bp_24_hour_summer <-
  bp_24_hour %>% 
  filter(season == "summer") %>% 
  mutate(across(where(is.double), ~as.numeric(scale(.)[,1])))

bp_24_hour_winter <-
  bp_24_hour %>% 
  filter(season == "winter") %>% 
  mutate(across(where(is.double), ~as.numeric(scale(.)[,1])))

#######################################################################################
# Top seasonal models without functional responses to burns.

# Individual- and population-level selection coefficients (conditional modes) for burns
# appear in Figure 3 and in Appendix S1 Table S5. 

# Fixed-effects coefficients for these models appear in Appendix S1 Table S4A.

# Note that 'map' and 'start' arguments tell glmmTMB to fix the variance of the random
# intercept at log(1e3) and freely estimate all nine random slopes in the models.
#######################################################################################

# TWO-WEEK MODELS
m_bp_2_week_summer_no_FR <-
  glmmTMB(case_ ~ tree  + TRI + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_2_week_summer, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))


m_bp_2_week_winter_no_FR <-
  glmmTMB(case_ ~ tree + I(tree^2) + TRI + I(TRI^2) + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_2_week_winter, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))


# 24 HOUR MODELS without functional responses

m_bp_24_hour_summer_no_FR <-
  glmmTMB(case_ ~ tree  + TRI + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_24_hour_summer, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))


m_bp_24_hour_winter_no_FR <-
  glmmTMB(case_ ~ tree + I(tree^2) + TRI + I(TRI^2) + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_24_hour_winter, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))

######################################################################################
# Top seasonal models including generalize functional response to burns.
# Predictions from these GFR models appear in Figure 4. 
######################################################################################

# TWO-WEEK MODELS
m_bp_2_week_summer <-
  glmmTMB(case_ ~ tree  + TRI + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            lc_burn:burn_avail_id + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_2_week_summer, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))

m_bp_2_week_winter <-
  glmmTMB(case_ ~ tree + I(tree^2) + TRI + I(TRI^2) + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            lc_burn:burn_avail_id + lc_burn:I(burn_avail_id^2) + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_2_week_winter, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))

saveRDS(m_bp_2_week_summer, "models/model_burn_perimeter_2_week_summer.rds")
saveRDS(m_bp_2_week_winter, "models/model_burn_perimeter_2_week_winter.rds")

# 24-HOUR MODELS
m_bp_24_hour_summer <-
  glmmTMB(case_ ~ tree  + TRI + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            lc_burn:burn_avail_id + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_24_hour_summer, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))

m_bp_24_hour_winter <-
  glmmTMB(case_ ~ tree + I(tree^2) + TRI + I(TRI^2) + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            lc_burn:burn_avail_id + lc_burn:I(burn_avail_id^2) + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_24_hour_winter, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))

saveRDS(m_bp_24_hour_summer, "models/model_burn_perimeter_24_hour_summer.rds")
saveRDS(m_bp_24_hour_winter, "models/model_burn_perimeter_24_hour_winter.rds")

#################################################################################
# Burn:ecotype interaction models, as in Appendix S1, Table S6
##################################################################################

# TWO-WEEK MODELS
m_bp_2_week_summer_ecotype <-
  glmmTMB(case_ ~ tree  + TRI + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            lc_burn:ecotype + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_2_week_summer, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))

m_bp_2_week_winter_ecotype <-
  glmmTMB(case_ ~ tree + I(tree^2) + TRI + I(TRI^2) + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            lc_burn:ecotype +  
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_2_week_winter, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))

# 24-HOUR MODELS
m_bp_24_hour_summer_ecotype <-
  glmmTMB(case_ ~ tree  + TRI + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            lc_burn:ecotype + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_24_hour_summer, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))

m_bp_24_hour_winter_ecotype <-
  glmmTMB(case_ ~ tree + I(tree^2) + TRI + I(TRI^2) + TPI + 
            lc_burn + lc_shrubs + lc_grass + lc_sparse + 
            lc_water + lc_other + lc_fen + log_step_length + 
            lc_burn:ecotype + 
            (1 | stratum) + (0 + tree | herd) + 
            (0 + TRI | herd) + (0 + lc_fen | herd) + (0 + lc_shrubs | herd) + 
            (0 + lc_other | herd) + (0 + lc_water | herd) + (0 + lc_grass | herd) + 
            (0 + lc_burn | herd/id), 
          data=bp_24_hour_winter, family=poisson,
          map=list(theta=factor(c(NA, 1:9))),
          start=list(theta=c(log(1e3), rep(0, 9))))

