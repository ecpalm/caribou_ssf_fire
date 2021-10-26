
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

# Script for running top models in within-burn RSA
# Code written by Eric Palm

############################################################################################

require(glmmTMB)
require(dplyr)

# Load data
wb <- readRDS("data/within_burn.rds") 

# Create seasonal datasets with centered and scaled continuous covariates
wb_summer <-
  wb %>% 
  filter(season == "summer") %>% 
  mutate(across(where(is.double), ~as.numeric(scale(.)[,1])))

wb_winter <-
  wb %>% 
  filter(season == "winter") %>% 
  mutate(across(where(is.double), ~as.numeric(scale(.)[,1])))

#######################################################################################
# Top seasonal models.

# Fixed-effects coefficients for these models appear in Appendix S1 Table S4C.

# Note that 'map' and 'start' arguments tell glmmTMB to fix the variance of the random
# intercept at log(1e3) and freely estimate all nine random slopes in the models.
#######################################################################################

# TWO-WEEK MODELS
m_wb_summer <-
  glmmTMB(case_ ~ d_edge + I(d_edge^2) + lichen + I(lichen^2) + severity:severity_01 + 
            TRI + I(TRI^2) + TPI + 
            lc_shrubs + lc_other + lc_grass + lc_sparse + lc_fen + 
            severity:severity_01:tsf + (1 | stratum) + 
            (0 + lc_shrubs | id) +
            (0 + severity:severity_01| id) + 
            (0 + lichen | id) + (0 + TRI | id), 
          data=wb_summer, family=poisson,
          map=list(theta=factor(c(NA, 1:4))),
          start=list(theta=c(log(1e3), rep(0, 4))))

m_wb_winter <-
  glmmTMB(case_ ~ d_edge + I(d_edge^2) + lichen + I(lichen^2) + severity:severity_01 + 
            TRI + I(TRI^2) + TPI + 
            lc_shrubs + lc_other + lc_grass + lc_sparse + lc_fen + 
            severity:severity_01:tsf + (1 | stratum) + (0 + d_edge|id) +
            (0 + lc_shrubs | id) + (0 + lc_fen |id) +
            (0 + severity:severity_01 | id) + 
            (0 + lichen | id) + (0 + TRI | id), 
          data=wb_winter, family=poisson,
          map=list(theta=factor(c(NA, 1:6))),
          start=list(theta=c(log(1e3), rep(0, 6))))


saveRDS(m_wb_summer, "models/model_within_burn_summer.rds")
saveRDS(m_wb_winter, "models/model_within_burn_winter.rds")
