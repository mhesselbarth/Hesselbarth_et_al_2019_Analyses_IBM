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
