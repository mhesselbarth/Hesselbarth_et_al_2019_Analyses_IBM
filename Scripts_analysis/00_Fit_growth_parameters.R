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

Rcpp::sourceCpp("Helper_functions/rcpp_calcl_ci.cpp", 
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
  dplyr::ungroup()

# plot growth vs dbh #
ggplot(beech_2013_top) +
  geom_point(aes(x = dbh_99, y = growth_full), pch = 1, size = 2) + 
  scale_x_continuous(name = "DBH 99 [cm]") +
  scale_y_continuous(name = "Mean anual growth [cm]") +
  theme_classic(base_size = 15)

# initialse function #
chapman <- growth_full ~ growth_assymp * growth_rate * growth_infl * 
  exp(-growth_rate * dbh_99) * ((1 - exp(-growth_rate * dbh_99)) ^ (growth_infl - 1))

# set starting functions adapted from Pommerening, A., Maleki, K., 2014. #
# Differences between competition kernels and traditional size-ratio based #
# competition indices used in forest ecology. For. Ecol. Manage. 331, 135-143. #
start_values <- list(growth_assymp = 75, 
                     growth_rate = 0.05, 
                     growth_infl = 3.5)

# fit model #
# fitted_fun <- nls(formula = chapman, data = beech_2013_top, start = start_values)
fitted_fun <- nlrq(formula = chapman, data = beech_2013_top, start = start_values,
                   trace = TRUE)

# predict values to plot function #
prediction_growth <- dplyr::bind_cols(dbh_99 = 1:100,
                                      growth_full = predict(object = fitted_fun,
                                                            data.frame(dbh_99 = 1:100)))

# plot result #
ggplot(beech_2013_top) +
  geom_point(aes(x = dbh_99, y = growth_full), pch = 1, size = 2) + 
  geom_line(data = prediction_growth,
            aes(x = dbh_99, y = growth_full), size = 1) +
  scale_x_continuous(name = "DBH 99 [cm]") +
  scale_y_continuous(name = "Mean anual growth [cm]") +
  theme_classic(base_size = 15)


# get summary of model fit #
# summary(fitted_fun)
broom::tidy(fitted_fun)

# # A tibble: 3 x 5
# term           estimate std.error statistic  p.value
# <chr>             <dbl>     <dbl>     <dbl>    <dbl>
# 1 growth_assymp 186.      27.2           6.81 2.59e-11
# 2 growth_rate     0.00686  0.000943      7.28 1.24e-12
# 3 growth_infl     1.54     0.0284       54.2  0.     

#### Fit ci parameters ####

# calculate sum distance to all neighbours
distance <- rcpp_calculate_ci(matrix = as.matrix(beech_2013[1:10, c("x", "y")]),
                              max_dist = 1e+9)

beech_2013_ci <- dplyr::mutate(beech_2013, distance = distance)

ggplot(beech_2013_ci) +
  geom_point(aes(x = dbh_99, y = growth_full), 
             pch = 1, size = 2, col = "#21908CFF") + 
  geom_point(data = beech_2013_top, 
             aes(x = dbh_99, y = growth_full), 
             pch = 1, size = 2, col = "#440154FF") + 
  geom_line(data = prediction_growth,
            aes(x = dbh_99, y = growth_full), size = 1) +
  scale_x_continuous(name = "DBH 99 [cm]") +
  scale_y_continuous(name = "Mean anual growth [cm]") +
  theme_classic(base_size = 15)

actual_growth <- growth_full ~ growth_pot * mod * ((1 - dbh_sum ^ beta * exp(-distance / dbh_sum ^ beta))) / dbh_99 ^ alpha + (((1 - dbh_sum ^ beta * exp(-distance / dbh_sum ^ beta))))

  
start_values <- list(growth_assymp = 75, 
                     growth_rate = 0.05, 
                     growth_infl = 3.5)

# actual_growth <- growth * growth_mod * )
