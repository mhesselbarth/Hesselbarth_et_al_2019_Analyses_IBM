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
NumericVector rcpp_calculate_ci(NumericMatrix matrix,
                                int max_dist) {
  
  // get number of rows
  int nrow = matrix.nrow();
  
  // initialise vector for ci value
  Rcpp::NumericVector distance(nrow, 0.0);
  
  // loop through all rows
  for(int i = 0; i < nrow - 1; i++){
    
    for(int j = i + 1; j < nrow; j++){
      
      // get distance between current point i and all points j
      const float dist_x = matrix(i, 0) - matrix(j, 0);
      const float dist_y = matrix(i, 1) - matrix(j, 1);
      
      const float distance_temp = std::sqrt(dist_x * dist_x + dist_y * dist_y);
      
      // distance above max_dist
      if(distance_temp > max_dist)
        continue; // nothing to do...
      
      // increase ci at point i and j
      distance[i] += distance_temp;
      distance[j] += distance_temp;
    }
  }
  return distance;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
mat <- matrix(runif(n = 200), ncol = 2)
mat_distance <- as.matrix(dist(mat, diag = TRUE, upper = TRUE))
mat_distance[mat_distance > 30] <- 0
as.numeric(apply(X = mat_distance, MARGIN = 1, FUN = sum))
rcpp_calculate_dist(mat, max_dist = 30)
*/
