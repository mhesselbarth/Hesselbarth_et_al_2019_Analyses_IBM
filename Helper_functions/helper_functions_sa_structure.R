#### SA_local_results_structure.R ####

calc_dbh_dist <- function(data, threshold, smaller = NULL, bigger = NULL, by = 1) {
  
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
    
    n_points <- nrow(data)
    
    message("> Progress: ", x, "/", length_input, " || Using ", n_points, " points.")
    
    dplyr::mutate(data, 
                  dbh_class = cut(dbh, breaks = seq(from = 0, 
                                                    to = max(dbh) + by, 
                                                    by = by), 
                                  labels = FALSE, include.lowest = TRUE, right = FALSE)) %>%
      dplyr::group_by(dbh_class) %>% 
      dplyr::summarise(n = dplyr::n(), 
                       n_rel = n / n_points) %>% 
      dplyr::mutate(dbh_class = dbh_class - 1) %>%
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
