###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Comparison spatial model structure ####

#### Import libraries and data ####

# load packages #
library(quantreg)
library(rabmp)
library(Rcpp)
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(spatstat)
library(tidyverse)

Rcpp::sourceCpp("Helper_functions/rcpp_calculate_actual.cpp", 
                embeddedR = FALSE)

# read paramters #
parameters_beech_default <- rabmp::read_parameters("Data/Input/parameters_beech.txt", return_list = TRUE)

# import data  #
pattern_2013 <- readr::read_rds("Data/Raw/pattern_2013_df.rds")

# filter data and calculate mean dbh growth #
beech_2013 <- dplyr::filter(pattern_2013, 
                            species == "beech", 
                            !is.na(dbh_99), 
                            !is.na(dbh_13),
                            type == "living", 
                            inside_fence == 0) %>% 
  dplyr::mutate(growth_mean = (growth_07 + growth_13) / 2,
                growth_full = (dbh_13 - dbh_99) / 14) %>% 
  dplyr::filter(growth_mean >= 0,
                growth_full >= 0)

#### Fit growth parameters ####

# compare growth calculated for full time period vs mean of 1999-2007 & 2007-2013 #
ggplot(beech_2013) +
  geom_point(aes(x = growth_mean, y = growth_full), pch = 1) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(name = "Growth mean 2007/2013", limits = c(-0.5, 1.5)) +
  scale_y_continuous(name = "Growth mean full time period", limits = c(-0.5, 1.5)) +
  coord_equal() + 
  theme_classic(base_size = 15)

cor(beech_2013$growth_mean, beech_2013$growth_full)

# classify into dbh classes and get only trees with highest growth #
breaks <- seq(from = 0, to = max(beech_2013$dbh_13) + 10, by = 10)

beech_2013_top <- dplyr::mutate(beech_2013,
                                dbh_class = cut(dbh_13, breaks = breaks)) %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::top_n(n = n() * 0.05, wt = growth_full) %>%
  # dplyr::top_n(n = 100, wt = growth_full) %>% 
  dplyr::pull(id)

# add top n classification to normal data
beech_2013 <- dplyr::mutate(beech_2013, 
                            top_n = dplyr::case_when(id %in% beech_2013_top ~ 1,
                                                     !id %in% beech_2013_top ~ 0), 
                            top_n = factor(top_n, levels = c(0, 1)))

# plot growth vs dbh #
ggplot(beech_2013) +
  geom_point(aes(x = dbh_99, y = growth_full, col = top_n), pch = 1, size = 2) + 
  scale_x_continuous(name = "DBH 99 [cm]") +
  scale_y_continuous(name = "Mean anual growth [cm]") +
  scale_color_viridis_d(name = "5% highest growth per class") +
  theme_classic(base_size = 15) + 
  theme(legend.position = "bottom")

# initialse function #
fun_potential <- function(dbh, assymp, rate, infl) {     
  
  growth <- assymp * rate * infl * exp(-rate * dbh) * 
    (1 - exp(-rate * dbh)) ^ (infl - 1)
  
  return(growth)
}

# set starting functions adapted from Pommerening, A., Maleki, K., 2014. #
# Differences between competition kernels and traditional size-ratio based #
# competition indices used in forest ecology. For. Ecol. Manage. 331, 135-143. #
start_values_potential <- list(assymp = 75, 
                               rate = 0.05, 
                               infl = 3.5)
          
# fit model #
fitted_fun_potential <- quantreg::nlrq(growth_full ~ 
                                         fun_potential(dbh_99, assymp, rate, infl), 
                         data = beech_2013, 
                         start = start_values_potential, 
                         tau = 0.995, 
                         trace = TRUE)

# predict values to plot function #
beech_2013 <- dplyr::mutate(beech_2013 , 
                            growth_pot = predict(object = fitted_fun_potential,
                                                 newdata = data.frame(dbh_99 = beech_2013$dbh_99)))

# plot result #
ggplot(beech_2013) +
  geom_point(aes(x = dbh_99, y = growth_full, col = top_n), pch = 1, size = 2) + 
  geom_line(aes(x = dbh_99, y = growth_pot), size = 1) +
  scale_x_continuous(name = "DBH 99 [cm]") +
  scale_y_continuous(name = "Mean anual growth [cm]") +
  scale_color_viridis_d(name = "5% highest growth per class") +
  theme_classic(base_size = 15) + 
  theme(legend.position = "bottom")

# get summary of model fit #
# summary(fitted_fun_potential)
broom::tidy(fitted_fun_potential)

# A tibble: 3 x 5
# term      estimate  std.error   statistic   p.value
# <chr>     <dbl>     <dbl>       <dbl>       <dbl>
# assymp    186.      28.8        6.45        1.17e-10
# rate      0.00732   0.00138     5.31        1.13e- 7
# infl      1.39      0.0512      27.2        0.     

#### Fit ci parameters ####

# initialse function #
fun_actual <- function(df, par) { 
  
  data_matrix <- as.matrix(df[, c("x", "y", "dbh_99", "growth_pot")])
  
  growth_modelled <- rcpp_calculate_actual(matrix = data_matrix, 
                                           modifier = par[1],
                                           alpha = par[2], 
                                           beta = par[3],
                                           max_dist = 30)
  
  difference <- sum(abs(df$growth_full - growth_modelled))
  
  return(difference)
}

# set starting functions adapted from Pommerening, A., Maleki, K., 2014. #
# Differences between competition kernels and traditional size-ratio based #
# competition indices used in forest ecology. For. Ecol. Manage. 331, 135-143. #
start_values_actual <- c(parameters_beech_default$growth_mod,
                         parameters_beech_default$ci_alpha, 
                         parameters_beech_default$ci_beta)

fitted_fun_actual <- optim(par = start_values_actual,
                           fn = fun_actual, 
                           df = beech_2013, 
                           method = "BFGS",
                           control = list(trace = TRUE, 
                                          maxit = 1000))

broom::tidy(fitted_fun_actual)
# A tibble: 3 x 2
# parameter   value
# <chr>       <dbl>
# parameter1  0.692
# parameter2  1.25 
# parameter3  0.336

fitted_fun_actual$value
# $value
# [1] 855.0206
