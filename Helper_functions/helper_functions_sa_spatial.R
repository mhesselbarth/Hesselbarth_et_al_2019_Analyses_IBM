#### Helper functions #### 

#### SA_local_results_spatial.R ####

calc_pcf <- function(data, window, r, 
                     fast = TRUE,
                     smaller = NULL, bigger = NULL, ...) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    if (is.null(smaller) && is.null(bigger)) {
      
      data <- dplyr::filter(data[[x]], i == max(i))
    }
    
    else{
      
      if (!is.null(smaller) && is.null(bigger)) { 
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh <= smaller)
      }
      
      else if (is.null(smaller) && !is.null(bigger)) {
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh >= bigger)
      }
      
      else if (!is.null(smaller) && !is.null(bigger)) {
        
        if (smaller <= bigger) {
        
          stop("Not possible to provide 'smaller' and 'bigger'.")
        }
        
        data <- dplyr::filter(data[[x]], i == max(i), 
                              dbh >= bigger, dbh <= smaller)
      }
    }
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
    
    data <- spatstat::ppp(x = data$x, y = data$y, 
                          window = window)
    
    if (!fast) {
      sf <- spatstat::pcf(data, r = r, ...)
    }
    
    else {
      sf <- shar::estimate_pcf_fast(data, r = r, ...)
    }
    
    tibble::add_column(tibble::as_tibble(sf), 
                       parameter = names_input[[x]], .before = 1)
    
  }, .id = "id") 
}

calc_nnd <- function(data, window, r, smaller = NULL, bigger = NULL, ...) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    if (is.null(smaller) && is.null(bigger)) {
      
      data <- dplyr::filter(data[[x]], i == max(i))
    }
    
    else{
      
      if (!is.null(smaller) && is.null(bigger)) { 
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh <= smaller)
      }
      
      else if (is.null(smaller) && !is.null(bigger)) {
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh >= bigger)
      }
      
      else if (!is.null(smaller) && !is.null(bigger)) {
        
        if (smaller <= bigger) {
          
          stop("Not possible to provide 'smaller' and 'bigger'.")
        }
        
        data <- dplyr::filter(data[[x]], i == max(i), 
                              dbh >= bigger, dbh <= smaller)
      }
    }
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
    
    data <- spatstat::ppp(x = data$x, y = data$y, 
                          window = window)
    
    sf <- spatstat::Gest(data, r = r, ...)
    
    tibble::add_column(tibble::as_tibble(sf), 
                       parameter = names_input[[x]], .before = 1)
    
  }, .id = "id") 
}

# calc_kmm <- function(data, window, length.out = 250, ...) {
#   
#   names_input <- names(data)
#   
#   length_input <- length(data)
#   
#   purrr::map_dfr(seq_along(data), function(x) {
#     
#     data <- dplyr::filter(data[[x]], i == max(i))
#     
#     message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
#     
#     data_ppp <- spatstat::ppp(x = data$x, y = data$y, 
#                               window = window)
#     
#     spatstat::marks(data_ppp) <- data$dbh
#     
#     r <- seq(from = 0,
#              to = spatstat::rmax.rule(W = data_ppp$window,
#                                       lambda = spatstat::intensity.ppp(data_ppp)),
#              length.out = length.out)
#     
#     sf <- spatstat::markcorr(data_ppp, ...)
#     
#     tibble::add_column(tibble::as_tibble(sf), 
#                        parameter = names_input[[x]], .before = 1)
#     
#   }, .id = "id") 
# }


calc_clark <- function(data, window, smaller = NULL, bigger = NULL, ...) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    if (is.null(smaller) && is.null(bigger)) {
      
      data <- dplyr::filter(data[[x]], i == max(i))
    }
    
    else{
      
      if (!is.null(smaller) && is.null(bigger)) { 
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh <= smaller)
      }
      
      else if (is.null(smaller) && !is.null(bigger)) {
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh >= bigger)
      }
      
      else if (!is.null(smaller) && !is.null(bigger)) {
        
        if (smaller <= bigger) {
          
          stop("Not possible to provide 'smaller' and 'bigger'.")
        }
        
        data <- dplyr::filter(data[[x]], i == max(i), 
                              dbh >= bigger, dbh <= smaller)
      }
    }
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
    
    data <- spatstat::ppp(x = data$x, y = data$y, 
                          window = window)
    
    index <- spatstat::clarkevans(data, ...)
    
    tibble::add_column(tibble::tibble(clark = index), 
                       parameter = names_input[[x]], .before = 1)
    
  }, .id = "id") 
}
