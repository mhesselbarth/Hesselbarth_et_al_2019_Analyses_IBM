###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Reconstruct input pattern ####

#### load packages ####
library(shar) # devtools::install_github("r-spatialecology/shar")
library(suppoRt) # devtools::install_github("mhesselbarth/suppoRt")
library(spatstat)
library(tidyverse)

source("Helper_functions/helper_functions_reconstruct_input.R")

#### import data ####
pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

# filter for all beech trees
beech_1999 <- spatstat::subset.ppp(pattern_1999, species == "beech")

#### reconstruct spatial pattern ####
beech_1999_rec <- reconstruction_helper(pattern = beech_1999,
                                        e_threshold = 0.001,
                                        max_runs = 20000, 
                                        select = "dbh_99")

# # get summary of results
print(beech_1999_rec, digits = 8)
shar::plot_energy(pattern = beech_1999_rec)
shar::plot_randomized_pattern(pattern = beech_1999_rec, ask = FALSE)

# # just to compare
# beech_1999_fit <- shar::fit_point_process(pattern = beech_1999, 
#                                           n_random = 1, 
#                                           process = "cluster")
# 
# print(beech_1999_fit)
# shar::calculate_energy(beech_1999_fit, return_mean = TRUE)
# shar::plot_energy(pattern = beech_1999_fit)
# shar::plot_randomized_pattern(pattern = beech_1999_fit, ask = FALSE)

# get only randomized pattern 
beech_1999_rec <- beech_1999_rec$randomized[[1]]

# update marks
spatstat::marks(beech_1999_rec) <- data.frame(species = "beech", type = "adult", 
                                              dbh = spatstat::marks(beech_1999_rec))


#### save results ####
suppoRt::save_rds(object = beech_1999_rec$randomized[[1]], 
                  filename = "beech_1999_rec.rds", 
                  path = "Data/Input/")
