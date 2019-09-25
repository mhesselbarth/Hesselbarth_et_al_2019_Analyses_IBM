###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Run biotic model ####

#### Import libraries and data ####

# load packages #
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(rabmp)
library(spatstat)
library(tidyverse)

parameters_beech_fitted <- rabmp::read_parameters("Data/Input/parameters_beech_fitted.txt", return_list = TRUE)

pattern_1999_recon <- readr::read_rds("Data/Input/beech_1999_rec_ppp.rds")

pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

#### Set SA parameters ####
repetitions <- 50 # 50

plot_area <- pattern_1999_recon$window
years <- rep(x = 50, times = repetitions) # 50
save_each <- 10
return_nested <- FALSE
verbose <- FALSE

#### Pre-processing of input data ####

data_reconstruction <- tibble::as_tibble(pattern_1999_recon) %>%
  dplyr::filter(species == "beech") %>%
  rabmp::prepare_data(x = "x", y = "y", species = "species", type = "type", dbh = "dbh")

data_real <- tibble::as_tibble(pattern_1999) %>%
  dplyr::select(x, y, species, dbh_99, type) %>% 
  dplyr::filter(species == "beech") %>%
  dplyr::mutate(type = "adult") %>% 
  rabmp::prepare_data(x = "x", y = "y", species = "species", type = "type", dbh = "dbh_99")

rm(pattern_1999_recon)

rm(pattern_1999)

#### Run model ####
# model_run_y100_r50_e5 <- purrr::map(years, function(x) {
#   rabmp::run_model(data = data,
#                    parameters = parameters_beech_default,
#                    plot_area = plot_area,
#                    years = x,
#                    save_each = save_each,
#                    return_nested = return_nested,
#                    verbose = TRUE)})

# reconstructed data #
model_run_y50_e10_r50_reco <- suppoRt::submit_to_cluster(rabmp::run_model,
                                                         years = years,
                                                         const = list(data = data_reconstruction,
                                                                      parameters = parameters_beech_fitted,
                                                                      plot_area = plot_area,
                                                                      save_each = save_each,
                                                                      return_nested = return_nested,
                                                                      verbose = verbose),
                                                         n_jobs = length(years),
                                                         template = list(job_name = "y50_e10_r50_reco",
                                                                         walltime = "12:00:00",
                                                                         queue = "medium", 
                                                                         mem_cpu = "4096", 
                                                                         log_file = "y50_e10_r50_reco.log"))

suppoRt::save_rds(object = model_run_y50_e10_r50_reco,
                  filename = "model_run_y50_e10_r50_reco.rds",
                  path = "Data/Output/")

rm(model_run_y50_e10_r50_reco)

# real world data #
model_run_y50_e10_r50_real <- suppoRt::submit_to_cluster(rabmp::run_model,
                                                         years = years,
                                                         const = list(data = data_real,
                                                                      parameters = parameters_beech_fitted,
                                                                      plot_area = plot_area,
                                                                      save_each = save_each,
                                                                      return_nested = return_nested,
                                                                      verbose = verbose),
                                                         n_jobs = length(years),
                                                         template = list(job_name = "y50_e10_r50_real",
                                                                         walltime = "12:00:00",
                                                                         queue = "medium", 
                                                                         mem_cpu = "4096", 
                                                                         log_file = "y50_e10_r50_real.log"))

suppoRt::save_rds(object = model_run_y50_e10_r50_real,
                  filename = "model_run_y50_e10_r50_real.rds",
                  path = "Data/Output/")

rm(model_run_y50_e10_r50_real)
