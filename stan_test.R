library(rstan)
library(Rcpp)
library(inline)

schools_dat <- list(J = 8, y = c(28,  8, -3,  7, -1,  1, 18, 12), sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
fit <- stan(file = '8schools.stan', data = schools_dat, iter = 1000, chains = 4)

print(fit, digits = 1)

strsplit(Sys.getenv('PATH'),';')

body <- ' NumericVector xx(x);
return wrap( std::accumulate( xx.begin(), xx.end(), 0.0));'

add <- inline::cxxfunction(signature(x = "numeric"), body, plugin = "Rcpp")

x <- 1
y <- 2
res <- add(c(x, y))
res

fx <- inline::cxxfunction(signature(x = "integer", y = "numeric"), '
  return  ScalarReal(  INTEGER(x)[0]  * 
  REAL(y)[0]  ) ;
')
fx(2L, 5)

devtools::find_rtools(TRUE)