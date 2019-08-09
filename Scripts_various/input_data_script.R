
# load packages
library(helpeR) # devtools::install_github("mhesselbarth/helpeR")
library(shar) # devtools::install_github("r-spatialecology/SHAR")
library(spatstat)
library(tidyverse)

# helper function for spatial and marks reconstruction
reconstruction_helper <- function(pattern, max_runs, select, species, type) {

  # reconstruction of spatial structure
  reconstructed <- shar::reconstruct_pattern_cluster(pattern = pattern,
                                                     max_runs = max_runs,
                                                     n_random = 1,
                                                     simplify =  TRUE,
                                                     return_input = FALSE,
                                                     verbose = FALSE)
  
  # reconstruction of DBH structure
  reconstructed <- shar::reconstruct_pattern_marks(pattern = reconstructed,
                                                   marked_pattern = spatstat::subset.ppp(pattern, select = select),
                                                   max_runs = max_runs,
                                                   n_random = 1,
                                                   simplify =  TRUE,
                                                   return_input = FALSE,
                                                   verbose = FALSE)

  # add species and DBH as marks to ppp
  spatstat::marks(reconstructed) <- data.frame(species = species, dbh = marks(reconstructed), type = type)

  return(reconstructed)
}

# import data
pattern_1999 <- readr::read_rds("Data/Raw/pattern_1999.rds")

table(pattern_1999$marks$Species)

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
max_runs <- 20000

# reconstruct all species and types 
pattern_1999_reconstructed <- runifpoint(n = 0, win = pattern_1999$window)

spatstat::marks(pattern_1999_reconstructed) <- data.frame(species = NA, dbh = NA, type = NA)

for (i in seq_along(seq_along(pattern_1999_split_flat))) {
  
  message("> Progress: ", i, "/", length(pattern_1999_split_flat), "\t\t\t", 
          appendLF = FALSE)
  
  name_split <- stringr::str_split(pattern_1999_split_flat_names[[i]], 
                                   pattern = "_", simplify = TRUE)
  
  current_reconstruction <- reconstruction_helper(pattern = pattern_1999_split_flat[[i]], 
                                                  species = name_split[1], type = name_split[2],
                                                  select = "DBH_99",
                                                  max_runs = max_runs)
  
  pattern_1999_reconstructed <- spatstat::superimpose(pattern_1999_reconstructed, 
                                                      current_reconstruction)
}

# save result
helpeR::save_rds(object = pattern_1999_reconstructed, 
                 filename = "pattern_1999_reconstructed.rds", 
                 path = "Data/Input", 
                 overwrite = FALSE)
