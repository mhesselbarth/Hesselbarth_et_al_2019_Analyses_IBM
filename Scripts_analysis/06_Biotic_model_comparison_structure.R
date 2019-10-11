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
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(spatstat)
library(tidyverse)

source("Helper_functions/helper_functions_comparison_structure.R")

pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")
pattern_2007 <- readr::read_rds("Data/Raw/pattern_2007_ppp.rds")
pattern_2013 <- readr::read_rds("Data/Raw/pattern_2013_ppp.rds")

model_run_y50_e5_r50 <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_real_b.rds")

names(model_run_y50_e5_r50) <- rep("Biotic model", 
                                   times = length(model_run_y50_e5_r50))

overwrite <- FALSE
base_size <- 10

#### Pre-processing data ####
df_1999 <- tibble::as_tibble(pattern_1999) %>% 
  dplyr::filter(species == "beech", 
                dbh_99 > 1)

df_2007 <- tibble::as_tibble(pattern_2007) %>% 
  dplyr::filter(species == "beech", 
                dbh_07 > 1, inside_fence == 0)

df_2013 <- tibble::as_tibble(pattern_2013) %>% 
  dplyr::filter(species == "beech", 
                dbh_13 > 1, inside_fence == 0)

#### Number trees alive ####
n_model <- calc_n_comp(data = model_run_y50_e5_r50)

n_2013 <- tibble::tibble(i = c(0, 8, 14), 
                         n = c(nrow(df_1999), nrow(df_2007), nrow(df_2013)))

#### Number ingrowth ####
n_ingrowth_model <- calc_n_ingrowth_comp(data = model_run_y50_e5_r50)

n_2013 <- tibble::tibble(i = c(8, 14), 
                         n = c(nrow(dplyr::filter(df_2007, is.na(dbh_99), 
                                                  !is.na(dbh_07), inside_fence == 0)),
                               nrow(dplyr::filter(df_2013, is.na(dbh_99), 
                                                  is.na(dbh_07), !is.na(dbh_13),
                                                  inside_fence == 0))))


#### Number trees dead #### BUG
n_dead_model <- calc_n_dead_comp(data = model_run_y50_e5_r50) %>% 
  dplyr::bind_rows() %>% 
  dplyr::group_by(i) %>% 
  dplyr::summarise(n_mean = mean(n), 
                   n_min = min(n), 
                   n_max = max(n), 
                   n_sd = sd(n))

n_dead_2013 <- tibble::tibble(i = c(8, 14), 
                              n = c(length(which(!df_1999$id %in% df_2007$id)), 
                                    length(which(!df_2007$id %in% df_2013$id))))

#### DBH distribution ####
# threshold <- 5
by <- 10

dbh_dist_model <- calc_dbh_dist_comp(data = model_run_y50_e5_r50, 
                                     by = by) %>% 
  dplyr::bind_rows() %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n_mean = mean(n), 
                   n_sd = sd(n),
                   n_rel_mean = mean(n_rel), 
                   n_rel_sd = sd(n_rel),) %>% 
  dplyr::mutate(data_type = "Biotic model")

dbh_dist_1999 <- dplyr::filter(df_1999, type != "dead", !is.na(dbh_99)) %>% 
  dplyr::mutate(dbh_class = cut(dbh_99, breaks = seq(from = 0, 
                                                     to = max(dbh_99) + by, 
                                                     by = by), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n_mean = dplyr::n(), 
                   n_rel_mean = n_mean / nrow(.)) %>% 
  dplyr::mutate(data_type = "Observed data 1999")

dbh_dist_2007 <- dplyr::filter(df_2007, 
                               type != "dead", !is.na(dbh_07), inside_fence == 0) %>% 
  dplyr::mutate(dbh_class = cut(dbh_07, breaks = seq(from = 0, 
                                                     to = max(dbh_07) + by, 
                                                     by = by), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n_mean = dplyr::n(), 
                   n_rel_mean = n_mean / nrow(.)) %>% 
  dplyr::mutate(data_type = "Observed data 2007")

dbh_dist_2013 <- dplyr::filter(df_2013, 
                               type != "dead", !is.na(dbh_13), inside_fence == 0) %>% 
  dplyr::mutate(dbh_class = cut(dbh_13, breaks = seq(from = 0, 
                                                     to = max(dbh_13) + by, 
                                                     by = by), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n_mean = dplyr::n(), 
                   n_rel_mean = n_mean / nrow(.)) %>% 
  dplyr::mutate(data_type = "Observed data 2013")

# dbh_dist_recon <- dplyr::filter(df_1999_reconstructed,
#                                 dbh > threshold, type != "dead") %>% 
#   dplyr::mutate(dbh_class = cut(dbh, breaks = seq(from = 0, 
#                                                      to = max(dbh) + by, 
#                                                      by = by), labels = FALSE)) %>%
#   dplyr::group_by(dbh_class) %>% 
#   dplyr::summarise(n = dplyr::n(), 
#                    n_rel = n / nrow(.)) %>% 
#   dplyr::mutate(dbh_class = dbh_class - 1) %>% 
#   tibble::add_column(parameter = "Reconstructed data", .before = 1)

dbh_dist_overall <- dplyr::bind_rows(dbh_dist_model, 
                                     dbh_dist_1999,
                                     dbh_dist_2007,
                                     dbh_dist_2013) %>% 
  dplyr::ungroup() %>% 
  tidyr::replace_na(replace = list(n_sd = 0, n_rel_sd = 0)) %>% 
  dplyr::mutate(data_type = factor(data_type, 
                                   levels = c("Biotic model", 
                                              "Observed data 1999", 
                                              "Observed data 2007",
                                              "Observed data 2013")), 
                dbh_class = factor(dbh_class, ordered = TRUE))

ggplot_biotic_dbh_dist <- ggplot(data = dbh_dist_overall) + 
  geom_bar(aes(x = dbh_class, y = n_rel_mean * 100, fill = data_type), 
           position = position_dodge(), stat = "identity") +
  scale_fill_viridis_d(name = "", option = "D") + 
  scale_x_discrete(name = "DBH class [cm]",
                   breaks = seq(from = 1, 
                                to = as.numeric(max(dbh_dist_overall$dbh_class)), 
                                by = 1), 
                   labels = paste0("<", seq(from = 1, 
                                            to = as.numeric(max(dbh_dist_overall$dbh_class)), 
                                            by = 1) * by)) +
  scale_y_continuous(name = "Relative frequency [%]",
                     breaks = seq(from = 0, to = 80, by = 10)) +
  
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_biotic_dbh_dist, 
                     filename = "ggplot_biotic_dbh_dist.png",
                     path = "Figures/Biotic_model/", 
                     dpi = 300, width = 15, height = 10, units = "cm", 
                     overwrite = overwrite)

#### DBH growth ####
by <- 20

min <- 0.1
low <- 0.25
median <- 0.5
high <- 0.75
max <- 0.9

dbh_growth_model <- calc_growth_comp(data = model_run_y50_e5_r50, by = by) %>% 
  dplyr::bind_rows() %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(inc_min = quantile(dbh_inc, probs = min), 
                   inc_low = quantile(dbh_inc, probs = low), 
                   inc_median = quantile(dbh_inc, probs = median),
                   inc_high = quantile(dbh_inc, probs = high),
                   inc_max = quantile(dbh_inc, probs = max)) %>% 
  dplyr::mutate(data_type = "Biotic model")

dbh_growth_2013 <- dplyr::filter(df_2013,  
                                 type != "dead", !is.na(dbh_13),!is.na(dbh_99),
                                 inside_fence == 0) %>% 
  dplyr::mutate(dbh_class = cut(dbh_99, breaks = seq(from = 0, 
                                                        to = max(dbh_99) + by, 
                                                        by = by), 
                                labels = FALSE)) %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(inc_min = quantile(growth_13, probs = min), 
                   inc_low = quantile(growth_13, probs = low), 
                   inc_median = quantile(growth_13, probs = median),
                   inc_high = quantile(growth_13, probs = high),
                   inc_max = quantile(growth_13, probs = max)) %>% 
  dplyr::mutate(data_type = "Observed data 2013")

dbh_growth_overall <- dplyr::bind_rows(dbh_growth_model, 
                                       dbh_growth_2013) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(data_type = factor(data_type, 
                                   levels = c("Biotic model", 
                                              "Observed data 2013")), 
                dbh_class = factor(dbh_class, ordered = TRUE))

ggplot_biotic_growth <- ggplot(data = dbh_growth_overall) +
  geom_boxplot(aes(x = dbh_class, fill = data_type,
                   ymin = inc_min, lower = inc_low, 
                   middle = inc_median,  
                   upper = inc_high, ymax = inc_max), 
               stat = "identity") + 
  # facet_wrap(~ data_type, scales = "free_y", ncol = 2, nrow = 2) +
  scale_fill_viridis_d(name = "") +
  scale_x_discrete(name = "DBH class [cm]",
                   breaks = seq(from = 1, 
                                to = as.numeric(max(dbh_growth_overall$dbh_class)), 
                                by = 1), 
                   labels = paste0("<", seq(from = 1, 
                                            to = as.numeric(max(dbh_growth_overall$dbh_class)), 
                                            by = 1) * by)) +
  scale_y_continuous(name = "Annual growth [cm]", 
                     breaks = seq(from = 0, to = max(dbh_growth_overall$inc_max),
                                   by = 0.2)) + 
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom")
  
suppoRt::save_ggplot(plot = ggplot_biotic_growth, 
                     filename = "ggplot_biotic_growth.png", 
                     path = "Figures/Biotic_model/", 
                     dpi = 300, width = 15, height = 12.5, units = "cm", 
                     overwrite = overwrite)
