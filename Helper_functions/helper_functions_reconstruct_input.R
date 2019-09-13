# helper function for spatial and marks reconstruction
reconstruction_helper <- function(pattern, max_runs = 10000, e_threshold = 0.01,
                                  select, verbose = TRUE) {
  
  # reconstruction of spatial structure
  reconstructed <- shar::reconstruct_pattern_cluster(pattern = pattern,
                                                     max_runs = max_runs,
                                                     e_threshold = e_threshold,
                                                     n_random = 1,
                                                     simplify =  TRUE,
                                                     return_input = FALSE,
                                                     verbose = verbose)
  
  
  marked_pattern <- spatstat::subset.ppp(pattern, select = select)
  
  # reconstruction of DBH structure
  reconstructed <- shar::reconstruct_pattern_marks(pattern = reconstructed,
                                                   marked_pattern = marked_pattern,
                                                   max_runs = max_runs,
                                                   e_threshold = e_threshold,
                                                   n_random = 1,
                                                   verbose = verbose)
  
  return(reconstructed)
}
