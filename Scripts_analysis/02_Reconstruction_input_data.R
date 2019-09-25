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

#### import data ####
pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999_ppp.rds")

# filter for all beech trees and DBH
beech_1999 <- spatstat::subset.ppp(pattern_1999, 
                                   species == "beech", select = dbh_99)

#### reconstruct spatial pattern ####
beech_1999_rec_rd_pat <- shar::reconstruct_pattern_cluster(pattern = beech_1999,
                                                    e_threshold = 0.001,
                                                    max_runs = 20000, 
                                                    n_random = 1,
                                                    simplify =  TRUE,
                                                    return_input = FALSE)

beech_1999_rec_rd_mar <- shar::reconstruct_pattern_marks(pattern = beech_1999_rec_rd_pat, 
                                                         marked_pattern = beech_1999,
                                                         max_runs = 20000,
                                                         e_threshold = 0.001,
                                                         n_random = 1) 

# # get summary of results
print(beech_1999_rec_rd_mar, digits = 8)
shar::plot_energy(pattern = beech_1999_rec_rd_mar)
shar::plot_randomized_pattern(pattern = beech_1999_rec_rd_mar, ask = FALSE)

suppoRt::save_rds(object = beech_1999_rec_rd_mar, 
                  filename = "beech_1999_rec_rd_mar.rds", 
                  path = "Data/Input/")

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
beech_1999_rec_ppp <- beech_1999_rec$randomized[[1]]

# update marks
spatstat::marks(beech_1999_rec_ppp) <- data.frame(species = "beech", type = "adult", 
                                                  dbh = spatstat::marks(beech_1999_rec_ppp))

#### save results ####
suppoRt::save_rds(object = beech_1999_rec_ppp, 
                  filename = "beech_1999_rec_ppp.rds", 
                  path = "Data/Input/")
