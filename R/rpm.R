# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Random Pivotal Method
#'
#' Selects samples with prescribed inclusion probabilities from a finite population.
#' This design has high entropy.
#' In each of the (at most) N steps, two undecided units are selected at random
#' to compete.
#'
#' If the inclusion probabilities \code{prob} sum to an integer n, the sample
#' size is fixed (n).
#'
#'
#' @param prob A vector of length N with inclusion probabilities.
#' @param eps A small value used to determine when an updated probability is
#' close enough to 0.0 or 1.0.
#'
#' @return A vector of selected indices in 1,2,...,N.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = rpm(prob);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#'
#' set.seed(12345);
#' prob = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9);
#' N = length(prob);
#' ep = rep(0L, N);
#' r = 10000L;
#' for (i in seq_len(r)) {
#'   s = rpm(prob);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
rpm = function(
  prob,
  eps = 1e-12
) {
  if (length(prob) == 1) {
    stop("'prob' must be a vector of probabilities");

    if (prob %% 1 != 0 || prob < 1)
      stop("'prob' must be a vector of probabilities or a single integer > 0");

    ## result = .rpm_int_cpp(prob, x, bucketSize, method);
  } else {
    if (length(prob) != dim(x)[2L])
      stop("the size of 'prob' and 'x' does not match");

    if (eps < 0.0 || 1e-4 < eps)
      stop("'eps' must be in [0.0, 1e-4]");

    result = .rpm_cpp(prob, eps);
  }

  return(result);
}
