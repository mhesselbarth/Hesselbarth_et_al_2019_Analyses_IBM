#### Import libraries and data ####

# load packages #
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(sensitivity)
library(spatstat)
library(tidyverse)

parameters_beech_default <- rabmp::read_parameters("Data/Input/parameters_beech.txt", return_list = TRUE)

pattern_1999_recon <- readr::read_rds("Data/Input/pattern_1999_reconstructed.rds")

source("Scripts/helper_functions.R")

#### Set SA parameters ####
repetitions <- 5 # 50
repetitions_hpc <- rep(x = 1, times = repetitions)

#### Create parameters ####
# increase parameters by 5% and 10%
parameters_beech_inc_5 <- change_parameters(x = parameters_beech_default, 
                                            change = 0.05) %>% 
  rep(each = repetitions)

parameters_beech_inc_10 <- change_parameters(x = parameters_beech_default, 
                                             change = 0.1) %>% 
  rep(each = repetitions)
# decrease parameters by 5% and 10%
parameters_beech_dec_5 <- change_parameters(x = parameters_beech_default, 
                                            change = -0.05) %>% 
  rep(each = repetitions)

parameters_beech_dec_10 <- change_parameters(x = parameters_beech_default, 
                                             change = -0.1) %>% 
  rep(each = repetitions)

parameters_beech_default <- list(parameters_beech_default) %>% 
  rep(each = repetitions)

#### Pre-processing of input data ####

data <- tibble::as_tibble(pattern_1999_recon) %>% 
  dplyr::filter(species == "Beech") %>% 
  dplyr::mutate(type = dplyr::case_when(type == "living" ~ "adult",
                                        type == "dead" ~ "dead"), 
                species = "beech") %>%
  rabmp::prepare_data(x = "x", y = "y", species = "species", type = "type", dbh = "dbh")

#### Set model arguments ####
plot_area <- pattern_1999_recon$window
years <- 50 # 25
seed_dispersal <- FALSE
save_each <- 5
verbose <- FALSE

#### Default parameters ####
sa_default <- suppoRt::submit_to_cluster(wrap_run_model, 
                                         parameters = parameters_beech_default,
                                         const = list(data = data, 
                                                      plot_area = plot_area, 
                                                      years = years, 
                                                      seed_dispersal = seed_dispersal,
                                                      save_each = save_each,
                                                      verbose = verbose),
                                         n_jobs = length(parameters_beech_default),
                                         template = list(job_name = "sa_default",
                                                         walltime = "00:30:00", 
                                                         queue = "medium"))

suppoRt::save_rds(object = sa_default,
                  filename = "sa_default.rds",
                  path = "Data/Output/")

# rm(sa_default)

#### Increased parameters #### 
sa_increased_5 <- suppoRt::submit_to_cluster(wrap_run_model, 
                                             parameters = parameters_beech_inc_5,
                                             const = list(data = data, 
                                                          plot_area = plot_area, 
                                                          years = years, 
                                                          seed_dispersal = seed_dispersal,
                                                          save_each = save_each,
                                                          verbose = verbose),
                                             n_jobs = length(parameters_beech_inc_5),
                                             template = list(job_name = "sa_inc_5",
                                                             walltime = "00:30:00", 
                                                             queue = "medium"))

suppoRt::save_rds(object = sa_increased_5,
                  filename = "sa_increased_5.rds",
                  path = "Data/Output/")

# rm(sa_increased_5)


sa_increased_10 <- suppoRt::submit_to_cluster(wrap_run_model, 
                                              parameters = parameters_beech_inc_10,
                                              const = list(data = data, 
                                                           plot_area = plot_area, 
                                                           years = years, 
                                                           seed_dispersal = seed_dispersal,
                                                           save_each = save_each,
                                                           verbose = verbose),
                                              n_jobs = length(parameters_beech_inc_10),
                                              template = list(job_name = "sa_inc_10",
                                                              walltime = "00:30:00", 
                                                              queue = "medium"))

suppoRt::save_rds(object = sa_increased_10,
                  filename = "sa_increased_10.rds",
                  path = "Data/Output/")

# rm(sa_increased_10)

#### Decreased parameters ####
sa_decreased_5 <- suppoRt::submit_to_cluster(wrap_run_model, 
                                             parameters = parameters_beech_dec_5,
                                             const = list(data = data, 
                                                          plot_area = plot_area, 
                                                          years = years, 
                                                          seed_dispersal = seed_dispersal,
                                                          save_each = save_each,
                                                          verbose = verbose),
                                             n_jobs = length(parameters_beech_dec_5),
                                             template = list(job_name = "sa_dec_5",
                                                             walltime = "00:30:00", 
                                                             queue = "medium"))

suppoRt::save_rds(object = sa_decreased_5,
                  filename = "sa_decreased_5.rds",
                  path = "Data/Output/")

# rm(sa_decreased_5)


sa_decreased_10 <- suppoRt::submit_to_cluster(wrap_run_model, 
                                              parameters = parameters_beech_dec_10,
                                              const = list(data = data, 
                                                           plot_area = plot_area, 
                                                           years = years, 
                                                           seed_dispersal = seed_dispersal,
                                                           save_each = save_each,
                                                           verbose = verbose),
                                              n_jobs = length(parameters_beech_dec_10),
                                              template = list(job_name = "sa_dec_10",
                                                              walltime = "00:30:00", 
                                                              queue = "medium"))

suppoRt::save_rds(object = sa_decreased_10,
                  filename = "sa_decreased_10.rds",
                  path = "Data/Output/")

# rm(sa_decreased_10)
