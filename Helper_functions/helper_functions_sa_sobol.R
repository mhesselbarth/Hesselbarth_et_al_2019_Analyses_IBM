###################################################
##    Author: Maximilian H.K. Hesselbarth        ##
##    Department of Ecosystem Modelling          ##
##    University of Goettingen                   ##
##    maximilian.hesselbarth@uni-goettingen.de   ##
##    www.github.com/mhesselbarth                ##
###################################################

#### Sobol indice ####
calc_sobol_indiv <- function(param, data, plot_area, years, save_each) {
  
  parameters <- list(ci_alpha = param[1], 
                     ci_beta = param[2],
                     ci_max_dist = parameters_fitted_biotic$ci_max_dist,
                     growth_assymp = parameters_fitted_biotic$growth_assymp, 
                     growth_infl = param[3],
                     growth_mod = parameters_fitted_biotic$growth_mod, 
                     growth_rate = parameters_fitted_biotic$growth_rate,
                     mort_dbh_early = param[4], 
                     mord_dbh_late = parameters_fitted_biotic$mort_dbh_late, 
                     mord_dinc = parameters_fitted_biotic$mort_dinc, 
                     mort_int_early = param[5], 
                     mort_int_late = parameters_fitted_biotic$mort_int_late, 
                     seed_str = parameters_fitted_biotic$seed_str, 
                     seed_success = parameters_fitted_biotic$seed_success,
                     seed_beta = parameters_fitted_biotic$seed_beta, 
                     seed_max_dist = parameters_fitted_biotic$seed_max_dist)
  
  output <- rabmp::run_model_biotic(data = data, 
                                    parameters = parameters, 
                                    plot_area = plot_area, 
                                    years = years,
                                    save_each = save_each, 
                                    verbose = FALSE)
  
  output <- dplyr::filter(output, i == max(i), type != "dead")
  
  nrow(output)
}

calc_sobol_pcf <- function(param) {
  
  parameters <- list(ci_alpha = param[1], 
                     ci_beta = param[2],
                     ci_max_dist = parameters_fitted_biotic$ci_max_dist,
                     growth_assymp = parameters_fitted_biotic$growth_assymp, 
                     growth_infl = param[3],
                     growth_mod = parameters_fitted_biotic$growth_mod, 
                     growth_rate = parameters_fitted_biotic$growth_rate,
                     mort_dbh_early = param[4], 
                     mord_dbh_late = parameters_fitted_biotic$mort_dbh_late, 
                     mord_dinc = parameters_fitted_biotic$mort_dinc, 
                     mort_int_early = parameters_fitted_biotic$mort_int_early, 
                     mort_int_late = param[5], 
                     seed_str = parameters_fitted_biotic$seed_str, 
                     seed_success = parameters_fitted_biotic$seed_success,
                     seed_beta = parameters_fitted_biotic$seed_beta, 
                     seed_max_dist = parameters_fitted_biotic$seed_max_dist)
  
  output <- rabmp::run_model_biotic(data = pattern_1999_dt, 
                                    parameters = parameters, 
                                    plot_area = plot_area, 
                                    years = years,
                                    save_each = save_each, 
                                    verbose = FALSE)
  
  output <- dplyr::filter(output, i == max(i), type != "dead")
  
  output <- spatstat::ppp(x = output$x, y = output$y, 
                          window = plot_area)
  
  output_pcf <- spatstat::pcf(output, correction = "Ripley", divisor = "d")
  
  MESS::auc(output_pcf$r, output_pcf$pcf, 
            from = min(output_pcf$r), to = max(output_pcf$r))
}
