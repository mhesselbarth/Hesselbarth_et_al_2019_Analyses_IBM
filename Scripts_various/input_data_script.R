
# load packages
library(UtilityFunctions) # devtools::install_github("mhesselbarth/UtilityFunctions")
library(SHAR) # devtools::install_github("r-spatialecology/SHAR")
library(spatstat)
library(tidyverse)

# helper function for spatial and marks reconstruction
reconstruction_helper <- function(pattern, select, species, type,
                                  n_random, max_runs,
                                  fitting, comp_fast,
                                  simplify, return_input) {

  # reconstruction of spatial structure
  x_reconstructed <- SHAR::reconstruct_pattern(pattern,
                                               n_random = n_random, max_runs = max_runs,
                                               fitting = fitting, comp_fast = comp_fast,
                                               simplify = simplify, return_input = return_input, 
                                               verbose = FALSE)
  
  # reconstruction of DBH structure
  x_reconstructed <- SHAR::reconstruct_marks(x_reconstructed,
                                             spatstat::subset.ppp(pattern, select = select),
                                             n_random = n_random, max_runs = max_runs,
                                             simplify = simplify, return_input = return_input,
                                             verbose = FALSE)

  # add species and DBH as marks to ppp
  spatstat::marks(x_reconstructed) <- data.frame(Species = species, DBH = marks(x_reconstructed), Type = type)

  return(x_reconstructed)
}


# import data
pattern_1999 <- readr::read_rds("C:/Users/Maximilian/ownCloud/08_Hainich/02_Modified_data/pattern_1999.rds")

table(pattern_1999$marks$Species)

# split data according to species (loop not possible because beech needs comp_fast = TRUE as only species)
# beech <- spatstat::subset.ppp(pattern_1999, Species == "Beech")
# ash <- spatstat::subset.ppp(pattern_1999, Species == "Ash")
# hornbeam <- spatstat::subset.ppp(pattern_1999, Species == "Hornbeam")
# sycamore <- spatstat::subset.ppp(pattern_1999, Species == "Sycamore")
# others <- spatstat::subset.ppp(pattern_1999, Species == "others")

# split data accordint to species and type
pattern_1999_split <- purrr::map(spatstat::split.ppp(pattern_1999, "Species"), 
                                 function(x) spatstat::split.ppp(x, "Type"))

pattern_1999_split_flat <- purrr::flatten(pattern_1999_split)

pattern_1999_split_flat_names <- c("Beech_dead", "Beech_living",
                                   "Ash_dead", "Ash_living",
                                   "Hornbeam_dead", "Hornbeam_living",
                                   "Sycamore_dead", "Sycamore_living",
                                   "others_dead", "others_living")

# set reconstruction parameters
n_random <- 1
max_runs <- 20000
fitting <- TRUE
comp_fast <- TRUE
simplify <- TRUE
return_input <- FALSE

# reconstruct all species and types 
pattern_1999_reconstructed <- runifpoint(n = 0, win = pattern_1999$window)
spatstat::marks(pattern_1999_reconstructed) <- data.frame(Species = NA, DBH = NA, Type = NA)

for(i in seq_along(seq_along(pattern_1999_split_flat))) {
  
  print(paste0("Progress: ", i, "/", length(pattern_1999_split_flat)))
  
  name_split <- stringr::str_split(pattern_1999_split_flat_names[[i]], 
                                   pattern = "_", simplify = TRUE)
  
  current_reconstruction <- reconstruction_helper(pattern = pattern_1999_split_flat[[i]], 
                                                  species = name_split[1], type = name_split[2],
                                                  select = "DBH_99",
                                                  n_random = n_random, max_runs = max_runs,
                                                  fitting = fitting, comp_fast = comp_fast,
                                                  simplify = simplify, return_input = return_input)
  
  pattern_1999_reconstructed <- spatstat::superimpose(pattern_1999_reconstructed, 
                                                      current_reconstruction)
}

# # reconstruct all species
# beech_reconstructed <- reconstruction_helper(pattern = beech, species = "Beech", select = "DBH_99",
#                                              n_random = n_random, max_runs = max_runs,
#                                              fitting = fitting, comp_fast = TRUE,
#                                              simplify = simplify, return_input = return_input)
# 
# ash_reconstructed <- reconstruction_helper(pattern = ash, species = "Ash", select = "DBH_99",
#                                            n_random = n_random, max_runs = max_runs,
#                                            fitting = fitting, comp_fast = FALSE,
#                                            simplify = simplify, return_input = return_input)
# 
# hornbeam_reconstructed <- reconstruction_helper(pattern = hornbeam, species = "Hornbeam", select = "DBH_99",
#                                                 n_random = n_random, max_runs = max_runs,
#                                                 fitting = fitting, comp_fast = FALSE,
#                                                 simplify = simplify, return_input = return_input)
# 
# sycamore_reconstructed <- reconstruction_helper(pattern = sycamore, species = "Sycamore", select = "DBH_99",
#                                                 n_random = n_random, max_runs = max_runs,
#                                                 fitting = fitting, comp_fast = FALSE,
#                                                 simplify = simplify, return_input = return_input)
# 
# others_reconstructed <- reconstruction_helper(pattern = others, species = "others", select = "DBH_99",
#                                               n_random = n_random, max_runs = max_runs,
#                                               fitting = fitting, comp_fast = FALSE,
#                                               simplify = simplify, return_input = return_input)
# 
# # superimpose all species to one pattern
# pattern_1999_reconstructed <- spatstat::superimpose(beech_reconstructed,
#                                                     ash_reconstructed,
#                                                     hornbeam_reconstructed,
#                                                     sycamore_reconstructed,
#                                                     others_reconstructed,
#                                                     W = pattern_1999$window)

# save result
UtilityFunctions::save_rds(pattern_1999_reconstructed, 
                           filename = "pattern_1999_reconstructed.rds", 
                           path = "C:/Users/Maximilian/ownCloud/03_Lehre/Msc_project_lena", 
                           overwrite = FALSE)
