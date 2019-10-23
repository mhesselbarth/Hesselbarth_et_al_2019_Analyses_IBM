###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Comparison model structure ####

#### Import libraries and data ####
source("Helper_functions/helper_functions_setup.R")
source("Helper_functions/helper_functions_comparison_structure.R")

pattern_2007 <- readr::read_rds("Data/Raw/pattern_2007_ppp.rds")
pattern_2013 <- readr::read_rds("Data/Raw/pattern_2013_ppp.rds")

model_run_y50_e5_r50_biotic <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_real_b.rds")
model_run_y50_e5_r50_abiotic <- readr::read_rds("Data/Output/model_runs/model_run_y50_e5_r50_real_a.rds")

#### Preprocess data ####

names(model_run_y50_e5_r50_biotic) <- rep("Biotic model", 
                                          times = length(model_run_y50_e5_r50_biotic))

names(model_run_y50_e5_r50_abiotic) <- rep("Abiotic model", 
                                           times = length(model_run_y50_e5_r50_abiotic))

df_2007 <- tibble::as_tibble(pattern_2007) %>% 
  dplyr::filter(species == "beech", 
                dbh_07 > 1, inside_fence == 0)

df_2013 <- tibble::as_tibble(pattern_2013) %>% 
  dplyr::filter(species == "beech", 
                dbh_13 > 1, inside_fence == 0)

#### DBH distribution ####
by_dist <- 10

dbh_dist_model_biotic <- calc_dbh_dist_comp(data = model_run_y50_e5_r50_biotic, 
                                            by = by_dist) %>% 
  dplyr::bind_rows() %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n_mean = mean(n), 
                   n_sd = sd(n),
                   n_rel_mean = mean(n_rel), 
                   n_rel_sd = sd(n_rel)) %>% 
  dplyr::mutate(data_type = "Biotic model")

dbh_dist_model_abiotic <- calc_dbh_dist_comp(data = model_run_y50_e5_r50_abiotic, 
                                             by = by_dist) %>% 
  dplyr::bind_rows() %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n_mean = mean(n), 
                   n_sd = sd(n),
                   n_rel_mean = mean(n_rel), 
                   n_rel_sd = sd(n_rel)) %>% 
  dplyr::mutate(data_type = "Abiotic model")

dbh_dist_2007 <- dplyr::filter(df_2007, 
                               type != "dead", !is.na(dbh_07), inside_fence == 0) %>% 
  dplyr::mutate(dbh_class = cut(dbh_07, breaks = seq(from = 0, 
                                                     to = max(dbh_07) + by_dist, 
                                                     by = by_dist), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n_mean = dplyr::n(), 
                   n_rel_mean = n_mean / nrow(.)) %>% 
  dplyr::mutate(data_type = "Field data 2007")

dbh_dist_2013 <- dplyr::filter(df_2013, 
                               type != "dead", !is.na(dbh_13), inside_fence == 0) %>% 
  dplyr::mutate(dbh_class = cut(dbh_13, breaks = seq(from = 0, 
                                                     to = max(dbh_13) + by_dist, 
                                                     by = by_dist), labels = FALSE)) %>%
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(n_mean = dplyr::n(), 
                   n_rel_mean = n_mean / nrow(.)) %>% 
  dplyr::mutate(data_type = "Field data 2013")

dbh_dist_overall <- dplyr::bind_rows(dbh_dist_model_biotic, 
                                     dbh_dist_model_abiotic,
                                     dbh_dist_2007,
                                     dbh_dist_2013) %>% 
  tidyr::replace_na(replace = list(n_sd = 0, n_rel_sd = 0)) %>% 
  dplyr::mutate(data_type = factor(data_type, 
                                   levels = c("Biotic model", 
                                              "Abiotic model", 
                                              "Field data 2007",
                                              "Field data 2013")), 
                dbh_class = factor(dbh_class, ordered = TRUE))

ggplot_dbh_dist <- ggplot(data = dbh_dist_overall) + 
  geom_bar(aes(x = dbh_class, y = n_rel_mean * 100, fill = data_type), 
           position = position_dodge(), stat = "identity") +
  scale_fill_viridis_d(name = "", option = "D") + 
  scale_x_discrete(name = "dbh class [cm]",
                   breaks = seq(from = 1, 
                                to = as.numeric(max(dbh_dist_overall$dbh_class)), 
                                by = 1), 
                   labels = paste0("<", seq(from = 1, 
                                            to = as.numeric(max(dbh_dist_overall$dbh_class)), 
                                            by = 1) * by_dist)) +
  scale_y_continuous(name = "Relative frequency [%]",
                     breaks = seq(from = 0, to = 60, by = 10),
                     limits = c(0, 60)) +
  theme_classic(base_size = base_size) + 
  theme(legend.position = "bottom", 
        legend.key.width = unit(0.5, units = "cm"))

suppoRt::save_ggplot(plot = ggplot_dbh_dist, 
                     filename = "ggplot_dbh_dist.png",
                     path = "Figures/", 
                     dpi = dpi, 
                     width = width_full, height = height_small, units = units, 
                     overwrite = overwrite)

#### Growth ####
by_growth <- 10

min <- 0.1
low <- 0.25
median <- 0.5
high <- 0.75
max <- 0.9

dbh_growth_model_biotic <- calc_growth_comp(data = model_run_y50_e5_r50_biotic, 
                                            by = by_growth) %>% 
  dplyr::bind_rows() %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(inc_min = quantile(dbh_inc, probs = min), 
                   inc_low = quantile(dbh_inc, probs = low), 
                   inc_median = quantile(dbh_inc, probs = median),
                   inc_high = quantile(dbh_inc, probs = high),
                   inc_max = quantile(dbh_inc, probs = max)) %>% 
  dplyr::mutate(data_type = "Biotic model")

dbh_growth_model_abiotic <- calc_growth_comp(data = model_run_y50_e5_r50_abiotic, 
                                             by = by_growth) %>% 
  dplyr::bind_rows() %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(inc_min = quantile(dbh_inc, probs = min), 
                   inc_low = quantile(dbh_inc, probs = low), 
                   inc_median = quantile(dbh_inc, probs = median),
                   inc_high = quantile(dbh_inc, probs = high),
                   inc_max = quantile(dbh_inc, probs = max)) %>% 
  dplyr::mutate(data_type = "Abiotic model")

dbh_growth_2013 <- dplyr::filter(df_2013,  
                                 type != "dead", !is.na(dbh_13),!is.na(dbh_99),
                                 inside_fence == 0) %>% 
  dplyr::mutate(dbh_class = cut(dbh_99, breaks = seq(from = 0, 
                                                     to = max(dbh_99) + by_growth, 
                                                     by = by_growth), 
                                labels = FALSE)) %>% 
  dplyr::group_by(dbh_class) %>% 
  dplyr::summarise(inc_min = quantile(growth_13, probs = min), 
                   inc_low = quantile(growth_13, probs = low), 
                   inc_median = quantile(growth_13, probs = median),
                   inc_high = quantile(growth_13, probs = high),
                   inc_max = quantile(growth_13, probs = max)) %>% 
  dplyr::mutate(data_type = "Field data 1999-2013")

dbh_growth_overall <- dplyr::bind_rows(dbh_growth_model_biotic,
                                       dbh_growth_model_abiotic,
                                       dbh_growth_2013) %>% 
  dplyr::mutate(data_type = factor(data_type, 
                                   levels = c("Biotic model",
                                              "Abiotic model",
                                              "Field data 1999-2013")), 
                dbh_class = factor(dbh_class, ordered = TRUE))

ggplot_growth <- ggplot(data = dbh_growth_overall) +
  geom_boxplot(aes(x = dbh_class, fill = data_type,
                   ymin = inc_min, lower = inc_low, 
                   middle = inc_median,  
                   upper = inc_high, ymax = inc_max), 
               stat = "identity") + 
  scale_fill_viridis_d(name = "") +
  scale_x_discrete(name = "dbh class [cm]",
                   breaks = seq(from = 1, 
                                to = as.numeric(max(dbh_growth_overall$dbh_class)), 
                                by = 1), 
                   labels = paste0("<", seq(from = 1, 
                                            to = as.numeric(max(dbh_growth_overall$dbh_class)), 
                                            by = 1) * by_growth)) +
  scale_y_continuous(name = "Annual growth [cm]", 
                     breaks = seq(from = 0, to = max(dbh_growth_overall$inc_max),
                                  by = 0.2)) + 
  theme_classic(base_size = base_size) +
  theme(legend.position = "bottom")

suppoRt::save_ggplot(plot = ggplot_growth, 
                     filename = "ggplot_growth.png", 
                     path = "Figures/", 
                     dpi = dpi, 
                     width = width_small, height = height_small, units = units, 
                     overwrite = overwrite)


