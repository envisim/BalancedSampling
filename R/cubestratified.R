# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Stratified Cube method with pooling of landing phases
#'
#' @inherit lpm1 description params return
#' @inherit cube params
#' @inheritSection cube Inclusion probabilities
#'
#' @param strata An integer vector of length N with stratum numbers.
#'
#' @references
#' Deville, J. C. and Tillé, Y. (2004).
#' Efficient balanced sampling: the cube method.
#' Biometrika, 91(4), 893-912.
#'
#' Chauvet, G. and Tillé, Y. (2006).
#' A fast algorithm for balanced sampling.
#' Computational Statistics, 21(1), 53-62.
#'
#' Chauvet, G. (2009).
#' Stratified balanced sampling.
#' Survey Methodology, 35, 115-119.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' strata = c(rep(1L, 100), rep(2L, 200), rep(3L, 300), rep(4L, 400));
#' s = cubestratified(prob, x, strata);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#'
#' set.seed(12345);
#' prob = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9);
#' N = length(prob);
#' x = matrix(runif(N * 2), ncol = 2);
#' strata = c(rep(1L, 1), rep(2L, 2), rep(3L, 3), rep(4L, 4));
#' ep = rep(0L, N);
#' r = 10000L;
#' for (i in seq_len(r)) {
#'   s = cubestratified(prob, x, strata);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
cubestratified = function(
  prob,
  x,
  strata,
  eps = 1e-12
) {
  if (!is.matrix(x)) {
    x = as.matrix(x);
  }

  if (length(strata) != dim(x)[1L])
    stop("the size of 'strata' and 'x' does not match");

  if (length(prob == 1))
    stop("'prob' must be a vector of probabilities");

  if (length(prob) != dim(x)[1L])
    stop("the size of 'prob' and 'x' does not match");

  if (eps < 0.0 || 1e-4 < eps)
    stop("'eps' must be in [0.0, 1e-4]");

  result = .cube_stratified_cpp(prob, x, strata, eps);

  return(result);
}
