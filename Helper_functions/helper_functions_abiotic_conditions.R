#### helper function abiotic conditions### 

add_padding <- function(raster, n = 1) {

  raster <- landscapemetrics::pad_raster(raster,
                                         pad_raster_value = NA, 
                                         pad_raster_cells = n)[[1]]
  
  names(raster) <- "value"
  
  cells_na <- raster::Which(is.na(raster),
                            cells = TRUE)
  
  neighbours <- raster::adjacent(x = raster,
                                 cells = cells_na,
                                 directions = 8)
  
  neighbours_value <- tibble::tibble(focal = neighbours[, 1],
                                     neighbour = neighbours[, 2],
                                     x = raster[neighbours[, 2]]) %>%
    dplyr::filter(!is.na(x)) %>%
    dplyr::group_by(neighbour, focal) %>%
    dplyr::summarise(x = mean(x, na.rm = TRUE))
  
  raster[neighbours_value$focal] <- neighbours_value$x
  
  return(raster)
}
