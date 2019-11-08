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
source("Helper_functions/helper_functions_setup.R")

Rcpp::sourceCpp("Helper_functions/rcpp_calculate_actual_biotic.cpp", 
                embeddedR = FALSE)

# read paramters #
parameters_default <- rabmp::read_parameters("Data/Input/parameters_default.txt", 
                                             sep = ";")

# import data  #
beech_2013 <- readr::read_rds("Data/Input/beech_2013_ppp.rds") %>% 
  tibble::as_tibble()

# filter data and calculate mean dbh growth #
beech_2013 <- dplyr::filter(beech_2013, 
                            !is.na(dbh_99), 
                            !is.na(dbh_13),
                            type == "living", 
                            inside_fence == 0) %>% 
  dplyr::mutate(growth_mean = (growth_07 + growth_13) / 2,
                growth_full = (dbh_13 - dbh_99) / 14) %>% 
  dplyr::filter(growth_mean >= 0,
                growth_full >= 0)

# compare growth calculated for full time period vs mean of 1999-2007 & 2007-2013 #
ggplot(beech_2013) +
  geom_point(aes(x = growth_mean, y = growth_full), pch = 1) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_x_continuous(name = "Growth mean 2007/2013", limits = c(-0.5, 1.5)) +
  scale_y_continuous(name = "Growth mean full time period", limits = c(-0.5, 1.5)) +
  coord_equal() + 
  theme_classic(base_size = base_size)

cor(beech_2013$growth_mean, beech_2013$growth_full)

#### Fit growth parameters ####

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
ggplot_fitting_dbh_growth <- ggplot(beech_2013) +
  geom_point(aes(x = dbh_99, y = growth_full, col = top_n), pch = 1, size = 2) + 
  scale_x_continuous(name = "dbh [cm]") +
  scale_y_continuous(name = "Mean annual growth [cm]", limits = c(0, 1.25)) +
  scale_color_manual(name = "5% highest growth per class", 
                     values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
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
start_values_potential <- list(assymp = 120, # max(beech_2013$dbh_99) 
                               rate = 0.005, # median(beech_2013$growth_full) / 14
                               infl = 1.5) # Try and error

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
ggplot_fitting_potential <- ggplot(beech_2013) +
  geom_point(aes(x = dbh_99, y = growth_full, col = top_n), pch = 1, size = 2) + 
  geom_line(aes(x = dbh_99, y = growth_pot), size = 1) +
  scale_x_continuous(name = "dbh [cm]") +
  scale_y_continuous(name = "Mean annual growth [cm]", limits = c(0, 1.25)) +
  scale_color_manual(name = "5% highest growth per class", 
                     values = c("#440154FF", "#21908CFF")) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom")

# get summary of model fit #
# summary(fitted_fun_potential)
broom::tidy(fitted_fun_potential)

# A tibble: 3 x 5
#   term      estimate  std.error   statistic   p.value
#   <chr>     <dbl>     <dbl>       <dbl>       <dbl>
#   assymp    205.      37.5        6.57        5.10e-11
#   rate      0.00649   0.00145     5.22        1.82e- 7
#   infl      1.35      0.0511      26.9        0.      

#### Fit ci parameters ####

# initialse function #
fun_actual <- function(df, par) { 
  
  data_matrix <- as.matrix(df[, c("x", "y", "dbh_99", "growth_pot")])
  
  growth_modelled <- rcpp_calculate_actual_biotic(matrix = data_matrix, 
                                                  alpha = par[1], 
                                                  beta = par[2], 
                                                  mod = 1,
                                                  max_dist = 30)
  
  difference <- sum(abs(df$growth_full - growth_modelled))
  
  return(difference)
}

# set starting functions adapted from Pommerening, A., Maleki, K., 2014. #
# Differences between competition kernels and traditional size-ratio based #
# competition indices used in forest ecology. For. Ecol. Manage. 331, 135-143. #
start_values_actual <- c(parameters_default$ci_alpha,
                         parameters_default$ci_beta)

# fit fun #
fitted_fun_actual <- optim(par = start_values_actual,
                           fn = fun_actual, 
                           df = dplyr::filter(beech_2013, 
                                              !id %in% beech_2013_top), 
                           method = "BFGS",
                           hessian = TRUE,
                           control = list(trace = TRUE, 
                                          maxit = 1000,
                                          REPORT = 1))

broom::tidy(fitted_fun_actual)
# A tibble: 3 x 2
# parameter   value
# <chr>       <dbl>
# parameter1  1.06 
# parameter2  0.435

fitted_fun_actual$value
# $value
# [1] 689.1816

standard_errors <- sqrt(abs(diag(solve(-fitted_fun_actual$hessian))))
# [1] 0.035125727 0.005524192

ci <- rabmp:::rcpp_calculate_ci(matrix = as.matrix(beech_2013[, c("x", "y", "dbh_99")]),
                                alpha = fitted_fun_actual$par[[1]],
                                beta = fitted_fun_actual$par[[2]],
                                max_dist = 30)

beech_2013 <- dplyr::mutate(beech_2013, ci = ci)

ggplot_fitting_actual <- ggplot(beech_2013) +
  geom_point(aes(x = dbh_99, y = growth_full, col = ci), pch = 1,  size = 2) + 
  geom_line(aes(x = dbh_99, y = growth_pot), size = 1) +
  scale_x_continuous(name = "DBH [cm]") +
  scale_y_continuous(name = "", limits = c(0, 1.25)) +
  scale_color_viridis_c(name = expression(c[i]^{trans}), option = "A") +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(1.5, "cm"))

#### Fit biotic seed dispersal ####
# Olesen, C.R., Madsen, P., 2008. The impact of roe deer (Capreolus capreolus),
# seedbed, light and seed fall on natural beech (Fagus sylvatica) regeneration.
# For. Ecol. Manag. 255, 3962–3972.

stand_1 <- 2.6 / 1167
stand_2 <- 3.7 / 652
stand_3 <- 2.6 / 307
default <- mean(c(stand_1, stand_2, stand_3))

#### Update parameters ####
parameters_fitted <- parameters_default

parameters_fitted$ci_alpha <- fitted_fun_actual$par[[1]]
parameters_fitted$ci_beta <- fitted_fun_actual$par[[2]]

parameters_fitted$growth_assymp <- broom::tidy(fitted_fun_potential)[[1, 2]]
parameters_fitted$growth_infl <- broom::tidy(fitted_fun_potential)[[3, 2]]
parameters_fitted$growth_mod <- 1
parameters_fitted$growth_rate <- broom::tidy(fitted_fun_potential)[[2, 2]]

parameters_fitted$seed_success <- default

write.table(parameters_fitted, row.names = FALSE, sep = ";")

#### Save plots #### 
ggplot_fitting_growth <- ggplot_fitting_potential + ggplot_fitting_actual + 
  patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")")

suppoRt::save_ggplot(plot = ggplot_fitting_growth, 
                     path = "Figures/",
                     filename = "ggplot_fitting_growth.png", 
                     dpi = dpi, units = units, 
                     height = height_full * 1/3, width = width_full, 
                     overwrite = overwrite)
