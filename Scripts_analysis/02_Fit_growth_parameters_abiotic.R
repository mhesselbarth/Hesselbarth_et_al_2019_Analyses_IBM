###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Fit model parameters ####

#### Import libraries and data ####
source("Helper_functions/helper_functions_setup.R")

# import helper function for fitting
Rcpp::sourceCpp("Helper_functions/rcpp_calculate_actual_abiotic.cpp",
                embeddedR = FALSE)

# read paramters #
parameters_default_abiotic <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt", 
                                                     sep = ";")

# set growth biotic parameter to 0
parameters_default_abiotic$growth_abiotic <- 0

# import data  #
abiotic_conditions_real <- readr::read_rds("Data/Input/abiotic_cond_real_fit.rds")

beech_2013_df <- readr::read_rds("Data/Input/beech_2013_ppp.rds") %>% 
  tibble::as_tibble()

#### Preprocessing of data ####

# filter data and calculate mean dbh growth #
beech_2013_df <- dplyr::filter(beech_2013_df, 
                               !is.na(dbh_99), 
                               !is.na(dbh_13),
                               type == "living", 
                               inside_fence == 0) %>% 
  dplyr::mutate(growth_mean = (growth_07 + growth_13) / 2,
                growth_full = (dbh_13 - dbh_99) / 14) %>% 
  dplyr::filter(growth_mean >= 0,
                growth_full >= 0)

# classify into dbh classes and get only trees with highest growth #
breaks <- seq(from = 0, to = max(beech_2013_df$dbh_13) + 10, by = 10)

beech_2013_df_top <- dplyr::mutate(beech_2013_df,
                                   dbh_class = cut(dbh_13, breaks = breaks)) %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::top_n(n = n() * 0.05, wt = growth_full) %>%
  dplyr::pull(id)

# initialse function potential growth #
fun_potential <- function(dbh, assymp, rate, infl) {     
  
  growth <- assymp * rate * infl * exp(-rate * dbh) * 
    (1 - exp(-rate * dbh)) ^ (infl - 1)
  
  return(growth)
}

# calculate potential growth #
beech_2013_df$growth_pot <- fun_potential(dbh = beech_2013_df$dbh_99, 
                                          assymp = parameters_default_abiotic$growth_assymp, 
                                          rate = parameters_default_abiotic$growth_rate, 
                                          infl = parameters_default_abiotic$growth_infl)

# set starting functions adapted from Pommerening, A., Maleki, K., 2014. #
# Differences between competition kernels and traditional size-ratio based #
# competition indices used in forest ecology. For. Ecol. Manage. 331, 135-143. #
start_values_actual <- c(ci_alpha = parameters_default_abiotic$ci_alpha, 
                         ci_beta = parameters_default_abiotic$ci_beta, 
                         growth_abiotic = parameters_default_abiotic$growth_abiotic)

#########################
####                 ####
#### Real world data ####
####                 ####
#########################

#### Get abiotic conditions ####
beech_2013_df$abiotic_real <- rabmp::extract_abiotic(data = data.table::data.table(x = beech_2013_df$x, 
                                                                                   y = beech_2013_df$y), 
                                                     abiotic = abiotic_conditions_real)[, 2]

#### Fit growth parameters ####

# initialse function #
fun_actual_real <- function(df, par) { 
  
  data_matrix <- as.matrix(df[, c("x", "y", "dbh_99", "growth_pot", "abiotic_real")])
  
  growth_modelled <- rcpp_calculate_actual_abiotic(matrix = data_matrix, 
                                                   alpha = par[1], 
                                                   beta = par[2],
                                                   mod = 1,
                                                   gamma = par[3],
                                                   max_dist = 30)
  
  difference <- sum(abs(df$growth_full - growth_modelled))
  
  return(difference)
}

# fit fun #
fitted_fun_actual_real <- optim(par = start_values_actual,
                                fn = fun_actual_real, 
                                df = dplyr::filter(beech_2013_df, 
                                                   !id %in% beech_2013_df_top), 
                                method = "BFGS",
                                hessian = TRUE,
                                control = list(trace = TRUE, 
                                               maxit = 1000,
                                               REPORT = 1))

broom::tidy(fitted_fun_actual_real)
# # A tibble: 3 x 3
# parameter      value std.error
# <chr>          <dbl>     <dbl>
# 1 ci_alpha       2.06    0.115  
# 2 ci_beta        0.286   0.00902
# 3 growth_abiotic 0.228   0.0294 

fitted_fun_actual_real$value
# [1] 669.8111

standard_errors <- sqrt(abs(diag(solve(-fitted_fun_actual_real$hessian))))
# ci_alpha      ci_beta       growth_abiotic 
# 0.115212436   0.009017715   0.029418124 

####################################
####                            ####
#### Seed dispersal & mortality ####
####                            ####
####################################

#### Fit abiotic seed dispersal ####
# Olesen, C.R., Madsen, P., 2008. The impact of roe deer (Capreolus capreolus),
# seedbed, light and seed fall on natural beech (Fagus sylvatica) regeneration.
# For. Ecol. Manag. 255, 3962–3972.

stand_1 <- 2.6 / 1167
stand_2 <- 3.7 / 652
stand_3 <- 2.6 / 307
default <- mean(c(stand_1, stand_2, stand_3))

# use CI from paper for abiotic
stand_1_hi <- 2.6 / 973 # 786.76
stand_2_hi <- 3.7 / 359 # 77.72
stand_3_hi <- 2.6 / 157 # 13 
high <- mean(c(stand_1_hi, stand_2_hi, stand_3_hi))

stand_1_lo <- 2.6 / 1361 # 1547.24
stand_2_lo <- 3.7 / 945 # 1226.28
stand_3_lo <- 2.6 / 457 # 601
low <- mean(c(stand_1_lo, stand_2_lo, stand_3_lo))

parameters_default_abiotic$seed_success_high <- high
parameters_default_abiotic$seed_success_low <- low

# Holzwarth, F., Kahl, A., Bauhus, J., Wirth, C., 2013. Many ways to die - 
# partitioning tree mortality dynamics in a near-natural mixed deciduous 
# forest. J. Ecol. 101, 220–230.

parameters_default_abiotic$mort_int_early_low <- 1.8 - ((1.8 - 1.2) / 1.96)
parameters_default_abiotic$mort_int_early_high <- 1.8 + ((2.5 - 1.8) / 1.96)

parameters_default_abiotic$mort_dbh_early_low <- -2.1 - ((-2.1 - -2.4) / 1.96)
parameters_default_abiotic$mort_dbh_early_high <- -2.1 + ((-1.9 - -2.1) / 1.96)

parameters_default_abiotic$mort_dinc_low <- -1.4 - ((-1.4 - -2.4) / 1.96)
parameters_default_abiotic$mort_dinc_high <- -1.4 + ((-0.45 - -1.4) / 1.96)

parameters_default_abiotic$mort_int_late_low <- -8.9 - ((-8.9 - -10.0) / 1.96)
parameters_default_abiotic$mort_int_late_high <- -8.9 + ((-7.8 - -8.9) / 1.96)

parameters_default_abiotic$mort_dbh_late_low <- 0.052 - ((0.052 - 0.033) / 1.96)
parameters_default_abiotic$mort_dbh_late_high <- 0.052 + ((0.070 - 0.052) / 1.96)

#### Update parameters ####
parameters_fitted_abiotic_real <- parameters_default_abiotic

parameters_fitted_abiotic_real$ci_alpha <- fitted_fun_actual_real$par[[1]]
parameters_fitted_abiotic_real$ci_beta <- fitted_fun_actual_real$par[[2]]
parameters_fitted_abiotic_real$growth_abiotic <- fitted_fun_actual_real$par[[3]]

#### Write parameters ####
write.table(parameters_fitted_abiotic_real, row.names = FALSE, sep = ";")
