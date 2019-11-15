###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Reconstruct input pattern ####

#### load packages ####
source("Helper_functions/helper_functions_setup.R")

#### import data ####
beech_1999_ppp <- readr::read_rds("Data/Input/beech_1999_ppp.rds") %>% 
  spatstat::subset.ppp(select = dbh_99)

#### reconstruct spatial pattern ####
beech_1999_ppp_rec <- shar::reconstruct_pattern_cluster(pattern = spatstat::unmark(beech_1999_ppp),
                                                        e_threshold = 0.001,
                                                        max_runs = 20000,
                                                        n_random = 1,
                                                        simplify =  TRUE,
                                                        return_input = FALSE)

beech_1999_ppp_rec <- shar::reconstruct_pattern_marks(pattern = beech_1999_ppp_rec, 
                                                      marked_pattern = beech_1999_ppp,
                                                      max_runs = 20000,
                                                      e_threshold = 0.001,
                                                      n_random = 1) 

# get summary of results
print(beech_1999_ppp_rec, digits = 8)
shar::plot_energy(pattern = beech_1999_ppp_rec)
shar::plot_randomized_pattern(pattern = beech_1999_ppp_rec, ask = FALSE)

# get only randomized pattern 
beech_1999_ppp_rec <- beech_1999_ppp_rec$randomized[[1]]

# update marks
spatstat::marks(beech_1999_ppp_rec) <- data.frame(species = "beech", type = "adult", 
                                                   dbh = spatstat::marks(beech_1999_ppp_rec))

#### save results ####
suppoRt::save_rds(object = beech_1999_ppp_rec, 
                  filename = "beech_1999_ppp_rec.rds", 
                  path = "Data/Input/", 
                  overwrite = overwrite)
