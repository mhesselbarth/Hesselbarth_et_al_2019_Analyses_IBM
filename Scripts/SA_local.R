#### Import libraries and data ####

# load packages #
library(suppoRt) # devtools::install_github("mhesselbarth/helpeR")
library(rabmp)
library(sensitivity)
library(spatstat)
library(tidyverse)

parameters_beech_default <- rabmp::read_parameters("Data/Input/parameters_beech.txt", return_list = TRUE)

pattern_1999_recon <- readr::read_rds("Data/Input/pattern_1999_reconstructed.rds")

#### Create parameters ####

# increase parameters by 5% and 10%
parameters_beech_inc_5 <- purrr::map(parameters_beech_default, function(x) x + x * 0.05) 
parameters_beech_inc_10 <- purrr::map(parameters_beech_default, function(x) x + x * 0.1) 

# decrease parameters by 5% and 10%
parameters_beech_dec_5 <- purrr::map(parameters_beech_default, function(x) x + x * 0.05) 
parameters_beech_dec_10 <- purrr::map(parameters_beech_default, function(x) x + x * 0.1) 


#### Pre-processing of input data ####

input_1999_beech <- tibble::as_tibble(pattern_1999_recon) %>% 
  dplyr::filter(species == "Beech") %>% 
  dplyr::mutate(type = dplyr::case_when(type == "living" ~ "adult",
                                        type == "dead" ~ "dead"), 
                species = "beech") %>%
  rabmp::prepare_data(x = "x", y = "y", species = "species", type = "type", dbh = "dbh")

#### Set SA parameters ####

repetitions <- 50
repetitions_hpc <- rep(x = 1, times = repetitions)
years <- 25

#### Model runs SA ####

# wrapper for HPC
wrapper_run_model <- function(x, ...) {
  
  rabmp::run_model(...)
}
  
# default #
# sa_default <- purrr::map(1:repetitions, function(x) run_model(data = input_1999_beech, 
#                                                               parameters = parameters_beech_default,
#                                                               plot_area = pattern_1999_recon$window, 
#                                                               years = years))

sa_default <- helpeR::submit_to_cluster(wrapper_run_model, 
                                        x = repetitions_hpc, 
                                        const = list(data = input_1999_beech, 
                                                     parameters = parameters_beech_default,
                                                     plot_area = pattern_1999_recon$window, 
                                                     years = years, 
                                                     verbose = FALSE),
                                        n_jobs = repetitions,
                                        template = list(job_name = "sa_default",
                                                        walltime = "01:30:00", 
                                                        queue = "fat",
                                                        mem_cpu = "20GB"))

suppoRt::save_rds(object = sa_default,
                 filename = "sa_default.rds",
                 path = "Data/Output/")

rm(sa_default)


