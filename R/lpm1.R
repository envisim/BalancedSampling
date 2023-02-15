# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Local Pivotal Method 1
#'
#' Selects spatially balanced samples with prescribed inclusion probabilities
#' from a finite population.
#'
#' Euclidean distance is used in the \code{x} space.
#'
#' If the inclusion probabilities \code{prob} sum to an integer n, the sample
#' size is fixed (n).
#'
#' The \code{type}s "kdtree" creates k-d-trees with terminal node bucket sizes
#' according to \code{bucketSize}.
#'
#' \itemize{
#'   \item{"kdtree0"} creates a k-d-tree using a median split on alternating variables.
#'   \item{"kdtree1"} creates a k-d-tree using a median split on the largest range.
#'   \item{"kdtree2"} creates a k-d-tree using a sliding-midpoint split.
#'   \item{"notree"} does a naive search for the nearest neighbour.
#' }
#'
#' @param prob A vector of length N with inclusion probabilities.
#' @param x An N by p matrix of (standardized) auxiliary variables.
#' @param type The method used in finding nearest neighbours.
#' Must be one of \code{"kdtree0"}, \code{"kdtree1"}, \code{"kdtree2"}, and
#' \code{"notree"}.
#' @param bucketSize The maximum size of the terminal nodes in the k-d-trees.
#' @param eps A small value used to determine when an updated probability is
#' close enough to 0.0 or 1.0.
#'
#' @return A vector of selected indices in 1,2,...,N.
#'
#' @references
#' Grafström, A., Lundström, N.L.P. & Schelin, L. (2012).
#' Spatially balanced sampling through the Pivotal method.
#' Biometrics 68(2), 514-520.
#'
#' Friedman, J. H., Bentley, J. L., & Finkel, R. A. (1977).
#' An algorithm for finding best matches in logarithmic expected time.
#' ACM Transactions on Mathematical Software (TOMS), 3(3), 209-226.
#'
#' Maneewongvatana, S., & Mount, D. M. (1999, December).
#' It’s okay to be skinny, if your friends are fat.
#' In Center for geometric computing 4th annual workshop on computational geometry (Vol. 2, pp. 1-8).
#'
#' Lisic, J. J., & Cruze, N. B. (2016, June).
#' Local pivotal methods for large surveys.
#' In Proceedings of the Fifth International Conference on Establishment Surveys.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = lpm1(prob, x);
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
#'   s = lpm1(prob, x);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
lpm1 = function(
  prob,
  x,
  type = "kdtree2",
  bucketSize = 50,
  eps = 1e-12
) {
  if (!is.matrix(x)) {
    x = t(as.matrix(x));
  } else {
    x = t(x);
  }

  if (type == "kdtree0") {
    method = 0;
  } else if (type == "kdtree1") {
    method = 1;
  } else if (type == "kdtree2") {
    method = 2;
  } else if (type == "notree") {
    method = 0;
    bucketSize = dim(x)[1L];
  } else {
    stop("'type' must be 'kdtree0', 'kdtree1', 'kdtree2', or 'notree'");
  }

  if (bucketSize %% 1 != 0 || bucketSize < 1)
    stop("'bucketSize' must be integer > 0");

  if (length(prob) == 1) {
    stop("'prob' must be a vector of probabilities");

    if (prob %% 1 != 0 || prob < 1)
      stop("'prob' must be a vector of probabilities or a single integer > 0");

    ## result = .lpm1_int_cpp(prob, x, bucketSize, method);
  } else {
    if (length(prob) != dim(x)[2L])
      stop("the size of 'prob' and 'x' does not match");

    if (eps < 0.0 || 1e-4 < eps)
      stop("'eps' must be in [0.0, 1e-4]");

    result = .lpm1_cpp(prob, x, bucketSize, method, eps);
  }

  return(result);
}
