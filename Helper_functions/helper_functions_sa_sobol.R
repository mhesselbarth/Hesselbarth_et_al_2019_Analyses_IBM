###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Sobol indice ####
calc_sobol_indiv <- function(x, data, parameters, plot_area, years, save_each) {
  
  parameters_temp <- list(ci_alpha = x[1], 
                          ci_beta = x[2],
                          ci_max_dist = parameters$ci_max_dist,
                          growth_assymp = parameters$growth_assymp, 
                          growth_infl = x[3],
                          growth_mod = parameters$growth_mod, 
                          growth_rate = parameters$growth_rate,
                          mort_dbh_early = x[4], 
                          mort_dbh_late = parameters$mort_dbh_late, 
                          mort_dinc = parameters$mort_dinc, 
                          mort_int_early = x[5], 
                          mort_int_late = x[6], 
                          seed_str = parameters$seed_str, 
                          seed_success = parameters$seed_success,
                          seed_eta = parameters$seed_eta, 
                          seed_max_dist = parameters$seed_max_dist)
  
  output <- rabmp::run_model_biotic(data = data, 
                                    parameters = parameters_temp, 
                                    plot_area = plot_area, 
                                    years = years,
                                    save_each = save_each, 
                                    verbose = FALSE)
  
  output <- dplyr::filter(output, i == max(i), type != "dead")
  
  output_adult <- dplyr::filter(output, type == "adult")
  
  output_sapling <- dplyr::filter(output, type == "sapling")
  
  list(adult = nrow(output_adult), sapling = nrow(output_sapling))
}

calc_sobol_pcf <- function(x, data, parameters, plot_area, years, save_each) {
  
  parameters_temp <- list(ci_alpha = x[1], 
                          ci_beta = x[2],
                          ci_max_dist = parameters$ci_max_dist,
                          growth_assymp = parameters$growth_assymp, 
                          growth_infl = x[3],
                          growth_mod = parameters$growth_mod, 
                          growth_rate = parameters$growth_rate,
                          mort_dbh_early = x[4], 
                          mort_dbh_late = parameters$mort_dbh_late, 
                          mort_dinc = parameters$mort_dinc, 
                          mort_int_early = x[5], 
                          mort_int_late = x[6], 
                          seed_str = parameters$seed_str, 
                          seed_success = parameters$seed_success,
                          seed_eta = parameters$seed_eta, 
                          seed_max_dist = parameters$seed_max_dist)
  
  output <- rabmp::run_model_biotic(data = data, 
                                    parameters = parameters_temp, 
                                    plot_area = plot_area, 
                                    years = years,
                                    save_each = save_each, 
                                    verbose = FALSE)
  
  output <- dplyr::filter(output, i == max(i), type != "dead")
  
  output_adult <- dplyr::filter(output, type == "adult")
  
  output_sapling <- dplyr::filter(output, type == "sapling")
  
  output_ppp_adult <- spatstat::ppp(x = output_adult$x, y = output_adult$y,
                                    window = plot_area)
  
  output_ppp_sapling <- spatstat::ppp(x = output_sapling$x, y = output_sapling$y,
                                      window = plot_area)
  
  output_pcf_adult <- spatstat::envelope(output_ppp_adult, 
                                         fun = "pcf", nsim = 199, verbose = FALSE,
                                         funargs = list(divisor = "d", correction = "good",
                                                        r = seq(from = 0, to = 75, 
                                                                length.out = 525)))
  
  output_pcf_sapling <- spatstat::envelope(output_ppp_sapling, 
                                           fun = "pcf", nsim = 199, verbose = FALSE,
                                           funargs = list(divisor = "d", correction = "good",
                                                          r = seq(from = 0, to = 75, 
                                                                  length.out = 525)))
  
  envl_summarised_adult <- onpoint::summarise_envelope(output_pcf_adult, 
                                                       seperated = FALSE, plot_result = FALSE)
  
  envl_summarised_sapling <- onpoint::summarise_envelope(output_pcf_sapling, 
                                                         seperated = FALSE, plot_result = FALSE)
  
  list(adult = envl_summarised_adult, sapling = envl_summarised_sapling)
}
