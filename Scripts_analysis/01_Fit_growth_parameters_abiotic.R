###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Fit model parameters ####

#### Import libraries and data ####

# load packages #
library(data.table)
library(rabmp)
library(Rcpp)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(spatstat)
library(tidyverse)

# import helper function for fitting
Rcpp::sourceCpp("Helper_functions/rcpp_calculate_actual_abiotic.cpp", 
                embeddedR = FALSE)

# read paramters #
parameters_default_abiotic <- rabmp::read_parameters("Data/Input/parameters_fitted_biotic.txt", 
                                                    sep = ";")

# set growth biotic parameter to 1
parameters_default_abiotic$growth_abiotic <- 1

# import data  #
abiotic_conditions <- readr::read_rds("Data/Input/abiotic_cond.rds")

pattern_2013_df <- readr::read_rds("Data/Raw/pattern_2013_df.rds")

pattern_1999_ppp <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

#### Preprocessing of data ####

# filter data and calculate mean dbh growth #
beech_2013_df <- dplyr::filter(pattern_2013_df, 
                               species == "beech", 
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

#### Get abiotic conditions ####
beech_2013_df$abiotic <- rabmp::extract_abiotic(data = data.table::data.table(x = beech_2013_df$x, 
                                                                              y = beech_2013_df$y), 
                                                abiotic = abiotic_conditions)

#### Fit growth parameters ####

# initialse function #
fun_actual <- function(df, par) { 
  
  data_matrix <- as.matrix(df[, c("x", "y", "dbh_99", "growth_pot", "abiotic")])
  
  growth_modelled <- rcpp_calculate_actual_abiotic(matrix = data_matrix, 
                                                   alpha = par[1], 
                                                   beta = par[2],
                                                   mod = 1,
                                                   gamma = par[3],
                                                   max_dist = 30)
  
  difference <- sum(abs(df$growth_full - growth_modelled))
  
  return(difference)
}

# set starting functions adapted from Pommerening, A., Maleki, K., 2014. #
# Differences between competition kernels and traditional size-ratio based #
# competition indices used in forest ecology. For. Ecol. Manage. 331, 135-143. #
start_values_actual <- c(parameters_default_abiotic$ci_alpha, 
                         parameters_default_abiotic$ci_beta, 
                         # parameters_default_abiotic$growth_mod,
                         parameters_default_abiotic$growth_abiotic)

# fit fun #
fitted_fun_actual <- optim(par = start_values_actual,
                           fn = fun_actual, 
                           df = dplyr::filter(beech_2013_df, 
                                              !id %in% beech_2013_df_top), 
                           method = "BFGS",
                           control = list(trace = TRUE, 
                                          maxit = 1000,
                                          REPORT = 1))

broom::tidy(fitted_fun_actual)
# A tibble: 3 x 2
# parameter   value
# <chr>       <dbl>
# parameter1  1.20 
# parameter2  0.359
# parameter3  4.10

fitted_fun_actual$value
# $value
# [1] 717.8145

# ci <- rabmp:::rcpp_calculate_ci(matrix = as.matrix(beech_2013_df[, c("x", "y", 
#                                                                      "dbh_99")]),
#                                 alpha = fitted_fun_actual$par[[1]],
#                                 beta = fitted_fun_actual$par[[2]],
#                                 max_dist = 30)
# 
# beech_2013_df <- dplyr::mutate(beech_2013_df, ci = ci)
# 
# ggplot_fitting_actual <- ggplot(beech_2013_df) +
#   geom_point(aes(x = dbh_99, y = growth_full, col = ci), pch = 1,  size = 2) + 
#   geom_line(aes(x = dbh_99, y = growth_pot), size = 1) +
#   scale_x_continuous(name = "DBH [cm]") +
#   scale_y_continuous(name = "Mean annual growth [cm]", limits = c(0, 1.25)) +
#   scale_color_viridis_c(name = "CI", option = "A") +
#   theme_classic(base_size = 15) + 
#   theme(legend.position = "bottom", 
#         legend.key.width = unit(1.5, "cm"))

#### Update parameters ####
parameters_fitted_abiotic <- parameters_default_abiotic

parameters_fitted_abiotic$ci_alpha <- fitted_fun_actual$par[[1]]
parameters_fitted_abiotic$ci_beta <- fitted_fun_actual$par[[2]]
parameters_fitted_abiotic$growth_abiotic <- fitted_fun_actual$par[[3]]

write.table(parameters_fitted_abiotic, row.names = FALSE, sep = ";")

#### Save plots #### 
# overwrite <- FALSE
# 
# suppoRt::save_ggplot(plot = ggplot_fitting_dbh_growth, 
#                      path = "Figures/Appendix",
#                      filename = "ggplot_fitting_dbh_growth.png", 
#                      dpi = 300, height = 10, width = 12.5, units = "cm", 
#                      overwrite = overwrite)
# 
# suppoRt::save_ggplot(plot = ggplot_fitting_potential, 
#                      path = "Figures/Appendix",
#                      filename = "ggplot_fitting_potential.png", 
#                      dpi = 300, height = 10, width = 12.5, units = "cm", 
#                      overwrite = overwrite)
# 
# suppoRt::save_ggplot(plot = ggplot_fitting_actual, 
#                      path = "Figures/Appendix",
#                      filename = "ggplot_fitting_actual.png", 
#                      dpi = 300, height = 10, width = 12.5, units = "cm", 
#                      overwrite = overwrite)


