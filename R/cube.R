# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' The Cube method
#'
#' @inherit lpm1 description params return
#'
#' @section Inclusion probabilities:
#' If the inclusion probabilities \code{prob} sum to an integer n, and the
#' inclusion probabilities are included as the first balancing variable,
#' the sample size is fixed (n).
#' prob, x eps, .cube_cpp, .cube_fast_cpp
#'
#' @param x An N by q matrix of balancing auxiliary variables.
#' @param fastFlight If \code{FALSE}, the flight phase will update all remaining
#' units in each step, which can be very time consuming for large populations.
#' Otherwise, it will utilize the fast flight method.
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

  if (length(prob == 1))
    stop("'prob' must be a vector of probabilities");

  if (length(prob) != dim(x)[1L])
    stop("the size of 'prob' and 'x' does not match");

  if (eps < 0.0 || 1e-4 < eps)
    stop("'eps' must be in [0.0, 1e-4]");

  if (fastFlight == FALSE) {
    result = .cube_cpp(prob, x, eps);
  } else {
    result = .cube_fast_cpp(prob, x, eps);
  }

  return(result);
}
