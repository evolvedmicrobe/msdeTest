#include "polarNormRand.h"

using namespace Rcpp;

double norm_keep = 0.0;

// Polar Form of Box-Muller Normal RNG See explanation at:
// http://www.design.caltech.edu/erik/Misc/Gaussian.html
// And it’s also explained on pages 494-496 of the A First Course in 
// Probability by Sheldon Ross, 7th Edition.
//[[Rcpp::export]]
double plnorm_rand() {
	double x1, x2, w, s;
	if(norm_keep != 0.0) { /* An exact test is intentional */
       s = norm_keep;
       norm_keep = 0.0;
       return s;
    } else {
		do {
	        x1 = 2.0 * unif_rand() - 1.0;
	        x2 = 2.0 * unif_rand() - 1.0;
	        w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
		w = sqrt( (-2.0 * LOG( w ) ) / w );
		norm_keep = x2 * w;
		return x1 * w;
	}
}

//' Simulate random standard normals using the polar form of the Box-Muller 
//' transformation.  Box-Muller is much faster than the Inversion method,
//' and this method is a faster version of that provided that the unif_rand
//' function is fast (otherwise performance is near BM).  Sample timings below.
//'
//'  > n = 5000000
//'  > message("R standard")
//'  R standard
//'  > RNGkind(kind="Mersenne-Twister")
//'  > RNGkind(normal.kind = "Inversion")
//'  > system.time(rnorm(n))
//'     user  system elapsed 
//'    0.444   0.004   0.450 
//'  > message("BM + Twister")
//'  BM + Twister
//'  > RNGkind(normal.kind = "Box-Muller")
//'  > system.time(rnorm(n))
//'     user  system elapsed 
//'    0.278   0.003   0.282 
//'  > message("Polar + Twister")
//'  Polar + Twister
//'  > system.time(msdeTest:::rplnorm(n))
//'     user  system elapsed 
//'    0.270   0.003   0.274 
//'  > message("Polar + Carry")
//'  Polar + Carry
//'  > RNGkind(kind="Marsaglia-Multicarry")
//'  > system.time(msdeTest:::rplnorm(n))
//'     user  system elapsed 
//'    0.234   0.003   0.238 
//' @param n Number of draws to take.
//' @return The size of the multinomial distributions
//' @export
// [[Rcpp::export]]
NumericVector rplnorm(int n) {
  NumericVector vals(n);
  for(int i = 0; i < n; i++) {
    vals[i] = plnorm_rand();
  }
  return vals;
}