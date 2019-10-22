###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Systematic find sigma ####

explore_sigma <- function(data,
                          sigma,
                          parameters,
                          probs,
                          probs_id,
                          plot_area,
                          years,
                          save_each) {
  
  # classify abiotic #
  
  # filter data using threshold sapling/adult #
  data_abiotic <- spatstat::subset.ppp(data, dbh_99 > 10 & species == "beech")
  
  # get intensity
  habitat <- spatstat::density.ppp(data_abiotic, at = "pixel",
                                   weights = data_abiotic$marks$dbh_99,
                                   dimyx = c(645, 609),
                                   kernel = "epanechnikov", sigma = sigma)
  
  # number of rows added at edges #
  n_rows <- 3
  
  # convert to raster and add padding #
  habitat <- tibble::as_tibble(habitat)
  
  habitat <- raster::rasterFromXYZ(habitat)
  
  habitat <- landscapemetrics::pad_raster(habitat,
                                          pad_raster_value = NA,
                                          pad_raster_cells = n_rows)[[1]]
  
  # add subsequently number of rows at edge #
  for (i in 1:n_rows) {
    
    # get all NA cells #
    cells_na <- raster::Which(is.na(habitat),
                              cells = TRUE)
    
    # get neigbors of NA cells #
    neighbours <- raster::adjacent(x = habitat,
                                   cells = cells_na,
                                   directions = 8)
    
    # get mean of all neighboring cells
    neighbours_value <- tibble::tibble(from = neighbours[, 1],
                                       to = neighbours[, 2],
                                       x_from = habitat[neighbours[, 1]],
                                       x_to = habitat[neighbours[, 2]])
    
    neighbours_value <- dplyr::group_by(neighbours_value, from)
    
    neighbours_value <- dplyr::summarise(neighbours_value,
                                         x = mean(x_to, na.rm = TRUE))
    
    neighbours_value <- dplyr::filter(neighbours_value, !is.na(x))
    
    # add values
    habitat[neighbours_value$from] <- neighbours_value$x
  }
  
  # scale value to -1 to 1 #
  habitat$scaled <- scales::rescale(raster::values(habitat),
                                    to = c(-1, 1), na.rm = TRUE)
  
  # set names #
  names(habitat) <- c("absolute", "scaled")
  
  # run model #
  data <- tibble::as_tibble(data)
  data <- dplyr::select(data, 
                        x, y, species, dbh_99, type)
  data <- dplyr::filter(data, species == "beech")
  data <- dplyr::mutate(data, type = "adult")
  data <- dplyr::select(data, -species)
  data <- rabmp::prepare_data(data,
                              x = "x", y = "y", type = "type", dbh = "dbh_99")
  
  result <-  rabmp::run_model_abiotic(data = data,
                                      parameters = parameters,
                                      abiotic = habitat$scaled,
                                      probs = probs[[probs_id]],
                                      plot_area = plot_area,
                                      years = years,
                                      save_each = save_each,
                                      verbose = FALSE)
  
  result$sigma <- sigma
  result$probs <- paste(probs[[probs_id]], collapse = "/")
  
  return(result)
}
