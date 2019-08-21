#### Helper functions #### 

#### SA_local.R ####
change_parameters <- function(x, change, return_list = TRUE) {
  
  parameters_beech_default_only <- unlist(x[which(x != 0)])
  
  parameters_others_default <- unlist(x[which(x == 0)])
  
  para <- purrr::map(seq_along(parameters_beech_default_only), function(i) {
    
    modification <- parameters_beech_default_only
    
    modification[i] <- modification[i] + modification[i] * change 
    
    modification <- c(modification, parameters_others_default)
    
    if (return_list) {
      
      modification <- as.list(modification)
    }
    
    return(modification)
  })
  
  names(para) <- names(parameters_beech_default_only)
  
  return(para)
}

#### SA_local_results_spatial.R ####

calc_pcf <- function(data, window, length.out = 250, ...) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    data <- dplyr::filter(data[[x]], i == max(i))
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
    
    data <- spatstat::ppp(x = data$x, y = data$y, 
                          window = window)
    
    r <- seq(from = 0,
             to = spatstat::rmax.rule(W = data$window,
                                      lambda = spatstat::intensity.ppp(data)),
             length.out = length.out)
    
    sf <- shar::estimate_pcf_fast(data, r = r, ...)
    
    tibble::add_column(tibble::as_tibble(sf), 
                       parameter = names_input[[x]], .before = 1)
    
  }, .id = "id") 
}

calc_nnd <- function(data, window, length.out = 250, ...) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    data <- dplyr::filter(data[[x]], i == max(i))
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
    
    data <- spatstat::ppp(x = data$x, y = data$y, 
                          window = window)
    
    r <- seq(from = 0,
             to = spatstat::rmax.rule(W = data$window,
                                      lambda = spatstat::intensity.ppp(data)),
             length.out = length.out)
    
    sf <- spatstat::Gest(data, ...)
    
    tibble::add_column(tibble::as_tibble(sf), 
                       parameter = names_input[[x]], .before = 1)
    
  }, .id = "id") 
}

calc_kmm <- function(data, window, length.out = 250, ...) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    data <- dplyr::filter(data[[x]], i == max(i))
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
    
    data_ppp <- spatstat::ppp(x = data$x, y = data$y, 
                          window = window)
    
    spatstat::marks(data_ppp) <- data$dbh
    
    r <- seq(from = 0,
             to = spatstat::rmax.rule(W = data_ppp$window,
                                      lambda = spatstat::intensity.ppp(data_ppp)),
             length.out = length.out)
    
    sf <- spatstat::markcorr(data_ppp, ...)
    
    tibble::add_column(tibble::as_tibble(sf), 
                       parameter = names_input[[x]], .before = 1)
    
  }, .id = "id") 
}
