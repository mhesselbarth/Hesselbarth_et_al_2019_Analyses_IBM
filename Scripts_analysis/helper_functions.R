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
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh < smaller)
      }
      
      else if (is.null(smaller) && !is.null(bigger)) {
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh > bigger)
      }
      
      else {
        
        stop("Not possible to provide 'smaller' and 'bigger'.")
      }
    }
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
    
    data <- spatstat::ppp(x = data$x, y = data$y, 
                          window = window)
    
    if (fast) {
      sf <- shar::estimate_pcf_fast(data, r = r, ...)
    }
    
    else{
      sf <- spatstat::pcf(data, r = r, ...)
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
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh < smaller)
      }
      
      else if (is.null(smaller) && !is.null(bigger)) {
        
        data <- dplyr::filter(data[[x]], i == max(i), dbh > bigger)
      }
      
      else {
        
        stop("Not possible to provide 'smaller' and 'bigger'.")
      }
    }
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
    
    data <- spatstat::ppp(x = data$x, y = data$y, 
                          window = window)
    
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


calc_clark <- function(data, window, ...) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    data <- dplyr::filter(data[[x]], i == max(i))
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data), " points.")
    
    data <- spatstat::ppp(x = data$x, y = data$y, 
                              window = window)
    
    index <- spatstat::clarkevans(data, ...)
    
    tibble::add_column(tibble::tibble(clark = index), 
                       parameter = names_input[[x]], .before = 1)
    
  }, .id = "id") 
}

#### SA_local_results_structure.R ####

calc_dbh_dist <- function(data, threshold, by = 1) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    data <- dplyr::filter(data[[x]], 
                          i == max(i), dbh > threshold, type != "dead")
    
    n_points <- nrow(data)
    
    message("> Progress: ", x, "/", length_input, " || Using ", n_points, " points.")
    
    dplyr::mutate(data, 
                  dbh_group = cut(dbh, breaks = seq(from = 0, 
                                                    to = max(dbh) + by, 
                                                    by = by), labels = FALSE)) %>%
      dplyr::group_by(dbh_group) %>% 
      dplyr::summarise(n = dplyr::n(), 
                       n_rel = n / n_points) %>% 
      dplyr::mutate(dbh_class = dbh_group - 1) %>%
      tibble::add_column(parameter = names_input[[x]], .before = 1)
  }, .id = "id")
}

calc_growth <- function(data, threshold, by = 1) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    years <- max(data[[x]]$i)
    
    data_start <- dplyr::filter(data[[x]], i == min(i), type != "dead")
    data_end <- dplyr::filter(data[[x]], i == max(i), type != "dead")
    
    data_merge <- dplyr::left_join(x = data_end[, c("id", "dbh")], 
                                   y = data_start[, c("id", "dbh")], 
                                   by = "id", suffix = c("_end", "_start"))
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data_merge), " points.")
    
    dplyr::mutate(data_merge, 
                  dbh_start = dplyr::case_when(is.na(dbh_start) ~ 0, 
                                   !is.na(dbh_start) ~ dbh_start),
                  dbh_inc = (dbh_end - dbh_start) / years) %>%
      tibble::add_column(parameter = names_input[[x]], .before = 1)
  }, .id = "id_map")
}

calc_died <- function(data) {
  
  names_input <- names(data)
  
  length_input <- length(data)
  
  purrr::map_dfr(seq_along(data), function(x) {
    
    data_dead <- dplyr::filter(data[[x]], type == "dead")
    
    data_start <- dplyr::filter(data[[x]], i == min(i))

    data_merge <- dplyr::anti_join(x = data_dead[, "id"], 
                                   y = data_start[, "id"], 
                                   by = "id", suffix = c("_dead", "_start"))
    
    message("> Progress: ", x, "/", length_input, " || Using ", nrow(data_merge), " points.")
    
    tibble::tibble(n_died = nrow(data_merge)) %>%
      tibble::add_column(parameter = names_input[[x]], .before = 1)
  }, .id = "id_map")
}
