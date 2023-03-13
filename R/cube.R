# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' The Cube method
#'
#' @description
#' Selects balanced samples with prescribed inclusion probabilities
#' from a finite population using the Cube Method.
#'
#' @details
#' If \code{prob} sum to an integer n, and \code{prob} is included as the first
#' balancing variable, a fixed sized sample (n) will be produced.
#'
#' @templateVar xbal x
#' @template sampling_template
#' @template x_template
#' @template probs_template
#'
#' @param fastFlight If `FALSE`, the flight phase will update all remaining
#' units in each step, which can be very time consuming for large populations.
#' Otherwise, it will use the fast flight method.
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
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = cube(prob, x);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#'
#' set.seed(12345);
#' prob = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9);
#' N = length(prob);
#' x = matrix(runif(N * 2), ncol = 2);
#' ep = rep(0L, N);
#' r = 10000L;
#' for (i in seq_len(r)) {
#'   s = cube(prob, x);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
cube = function(
  prob,
  x,
  fastFlight = TRUE,
  eps = 1e-12
) {
  if (!is.matrix(x)) {
    x = as.matrix(x);
  }

  N = dim(x)[1L];
  .eps_check(eps);
  prob = .prob_check(prob, N);

  if (fastFlight == FALSE) {
    result = .cube_cpp(prob, x, eps);
  } else {
    result = .cube_fast_cpp(prob, x, eps);
  }

  return(result);
}
