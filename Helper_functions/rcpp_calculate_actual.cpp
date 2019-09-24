#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector rcpp_calculate_actual(NumericMatrix matrix,
                                    double alpha, double beta,
                                    int max_dist) {
  
  // get number of rows
  const int nrow = matrix.nrow();
  
  // initialise matrix for ci value and dbh
  NumericVector result(nrow, 0.0);
  
  // initialise vector for ci value
  Rcpp::NumericVector ci(nrow, 0.0);
  
  // loop through all rows
  for(int i = 0; i < nrow - 1; i++){
    
    // initialse float for dbh
    const float dbh_i = matrix(i, 2);
    const float pot_i = matrix(i, 3);
    
    for(int j = i + 1; j < nrow; j++){
      
      const float dbh_j = matrix(j, 2);
      
      // get distance between current point i and all points j
      const float dist_x = matrix(i, 0) - matrix(j, 0);
      const float dist_y = matrix(i, 1) - matrix(j, 1);
      
      const float distance = std::sqrt(dist_x * dist_x + dist_y * dist_y);
      
      // distance above max_dist
      if(distance > max_dist)
        continue; // nothing to do...
      
      // calculate ci of current i and j
      const float ci_temp_i = std::pow(dbh_j, alpha) * std::exp(-distance / std::pow(dbh_j, beta));
      
      const float ci_temp_j = std::pow(dbh_i, alpha) * std::exp(-distance / std::pow(dbh_i, beta));
      
      // increase ci at point i and j
      ci[i] += ci_temp_i;
      ci[j] += ci_temp_j;
    }
    
    // normalize between 0 - 1 and scale with own dbh
    ci[i] = ci[i] / (std::pow(dbh_i, alpha) + ci[i]);
    
    // calculate actual growth
    // result[i] = pot_i * modifier * (1 - ci[i]);
    result[i] = pot_i * (1 - ci[i]);
  }
  
  // normalize last ci
  ci[nrow - 1] = ci[nrow - 1] / (std::pow(matrix(nrow - 1, 2), alpha) + ci[nrow - 1]);
  
  // calculate actual growth
  // result[nrow - 1] = matrix(nrow - 1, 3) * modifier * (1 - ci[nrow - 1]);
  result[nrow - 1] = matrix(nrow - 1, 3) * (1 - ci[nrow - 1]);
  
  return result;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# Add example
*/
