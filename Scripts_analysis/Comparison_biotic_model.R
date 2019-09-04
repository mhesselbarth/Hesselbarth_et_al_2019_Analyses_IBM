#### Import libraries and data ####

# load packages #
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(spatstat)
library(tidyverse)

source("Scripts_analysis/helper_functions.R")

pattern_1999 <- readr::read_rds("Data/Input/pattern_1999.rds")
pattern_2007 <- readr::read_rds("Data/Input/pattern_2007.rds")
pattern_2013 <- readr::read_rds("Data/Input/pattern_2013.rds")

# pattern_1999_reconstructed <- readr::read_rds("Data/Input/pattern_1999_reconstructed.rds")

model_run_y100_r50_e100 <- readr::read_rds("Data/Output/model_run_y100_r50_e100.rds")

names(model_run_y100_r50_e100) <- rep("Biotic model", times = length(model_run_y100_r50_e100))

# source("Scripts_analysis/helper_functions.R")

overwrite <- FALSE

#### Pre-processing data ####
df_1999 <- tibble::as_tibble(pattern_1999)
df_2007 <- tibble::as_tibble(pattern_2007)
df_2013 <- tibble::as_tibble(pattern_2013)
# df_1999_reconstructed <- tibble::as_tibble(pattern_1999_reconstructed)

#### DBH distribution ####
threshold <- 5
by <- 1

dbh_dist_model <- calc_dbh_dist(data = model_run_y100_r50_e100, 
                                threshold = threshold, by = by) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(n = mean(n), 
                   n_rel = min(n_rel))

dbh_dist_1999 <- dplyr::filter(df_1999,
                               DBH_99 > threshold, Type != "dead") %>% 
  dplyr::mutate(dbh_class = cut(DBH_99, breaks = seq(from = 0, 
                                                     to = max(DBH_99) + by, 
                                                     by = by), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n = dplyr::n(), 
                   n_rel = n / nrow(.)) %>% 
  dplyr::mutate(dbh_class = dbh_class - 1) %>% 
  tibble::add_column(parameter = "Observed data 1999", .before = 1)

dbh_dist_2007 <- dplyr::filter(df_2007,
                               DBH_07 > threshold, Type != "dead") %>% 
  dplyr::mutate(dbh_class = cut(DBH_07, breaks = seq(from = 0, 
                                                     to = max(DBH_07) + by, 
                                                     by = by), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n = dplyr::n(), 
                   n_rel = n / nrow(.)) %>% 
  dplyr::mutate(dbh_class = dbh_class - 1) %>% 
  tibble::add_column(parameter = "Observed data 2007", .before = 1)

dbh_dist_2013 <- dplyr::filter(df_2013,
                               DBH_13 > threshold, Type != "dead") %>% 
  dplyr::mutate(dbh_class = cut(DBH_13, breaks = seq(from = 0, 
                                                     to = max(DBH_13) + by, 
                                                     by = by), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n = dplyr::n(), 
                   n_rel = n / nrow(.)) %>% 
  dplyr::mutate(dbh_class = dbh_class - 1) %>% 
  tibble::add_column(parameter = "Observed data 2013", .before = 1)

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
  dplyr::mutate(parameter = forcats::as_factor(parameter))

ggplot_dbh_dist <- ggplot(data = dbh_dist_overall) + 
  geom_bar(aes(x = dbh_class, y = n_rel, fill = parameter), 
           position = "dodge2", stat = "identity") +
  scale_fill_viridis_d(name = "Data type", option = "D") +
  scale_x_continuous(name = "DBH class [cm]",
                     breaks = seq(from = threshold,
                                  to = 115,
                                  by = 10),
                     labels = paste0(">", seq(from = threshold,
                                              to = 115,
                                              by = 10))) +
  coord_cartesian(xlim = c(threshold, 115)) +
  scale_y_continuous(name = "Relative frequency") +
  theme_classic(base_size = 15)

suppoRt::save_ggplot(plot = ggplot_dbh_dist, 
                     filename = "ggplot_dbh_dist.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)

#### DBH growth ####
by <- 20

dbh_growth_model <- calc_growth(data = model_run_y100_r50_e100) %>% 
  dplyr::mutate(dbh_class = cut(dbh_start, breaks = seq(from = 0, 
                                                        to = max(dbh_start) + by, 
                                                        by = by), labels = FALSE, 
                                include.lowest = TRUE)) %>% 
  dplyr::group_by(parameter, dbh_class) %>% 
  dplyr::summarise(min = min(dbh_inc), 
                   lower = quantile(dbh_inc, 0.25), 
                   median = median(dbh_inc), 
                   higher =  quantile(dbh_inc, 0.75),
                   max = max(dbh_inc))

dbh_growth_2013 <- dplyr::filter(df_2013, Type != "dead") %>%
  dplyr::mutate(DBH_99 = dplyr::case_when(is.na(DBH_99) ~ 0, 
                                          !is.na(DBH_99) ~ DBH_99),
                dbh_inc = (DBH_13 - DBH_99) / 14, 
                dbh_class = cut(DBH_99, breaks = seq(from = 0, 
                                                     to = max(DBH_99, na.rm = TRUE) + by,
                                                     by = by), labels = FALSE, 
                                include.lowest = TRUE)) %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(min = min(dbh_inc), 
                   lower = quantile(dbh_inc, 0.25), 
                   median = median(dbh_inc), 
                   higher =  quantile(dbh_inc, 0.75),
                   max = max(dbh_inc)) %>% 
  tibble::add_column(parameter = "Observed data 2013", .before = 1)

dbh_growth_overall <- dplyr::bind_rows(dbh_growth_model, 
                                       dbh_growth_2013) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(parameter = forcats::as_factor(parameter), 
                dbh_class = forcats::as_factor(dbh_class))

ggplot_growth <- ggplot(data = dbh_growth_overall) +
  geom_boxplot(aes(x = dbh_class, 
                   ymin = min, lower = lower, 
                   middle = median, 
                   upper = higher, ymax = max, 
                   fill = parameter), 
               stat = "identity") + 
  # facet_wrap(~ parameter) +
  # scale_x_continuous(labels = paste0(">", seq(from = 1,
  #                                             to = 115,
  #                                             by = 10))) +
  theme_classic(base_size = 15)
  
  
suppoRt::save_ggplot(plot = sa_ggplot_growth, 
                     filename = "sa_ggplot_growth.png", path = "Figures/", 
                     dpi = 300, width = 30, height = 15, units = "cm", 
                     overwrite = overwrite)

