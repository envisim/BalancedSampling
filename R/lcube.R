# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' The Local Cube method
#'
#' @description
#' Selects doubly balanced samples with prescribed inclusion probabilities
#' from a finite population.
#'
#' @details
#' Euclidean distance is used in the \code{Xspread} space.
#'
#' @inherit lpm1 params return
#' @inheritSection cube Inclusion probabilities
#' @inheritSection lpm1 k-d-trees
#'
#' @param Xspread An N by p matrix of (standardized) auxiliary variables.
#' @param Xbal An N by q matrix of balancing auxiliary variables.
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
#' Grafström, A. and Tillé, Y. (2013).
#' Doubly balanced spatial sampling with spreading and restitution of auxiliary totals.
#' Environmetrics, 24(2), 120-131
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = lcube(prob, x, prob);
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
#'   s = lcube(prob, x, prob);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
lcube = function(
  prob,
  Xspread,
  Xbal,
  type = "kdtree2",
  bucketSize = 50,
  eps = 1e-12
) {
  if (!is.matrix(Xbal)) {
    Xbal = as.matrix(Xbal);
  }

  if (!is.matrix(Xspread)) {
    Xspread = t(as.matrix(Xspread));
  } else {
    Xspread = t(Xspread);
  }

  if (type == "kdtree0") {
    method = 0;
  } else if (type == "kdtree1") {
    method = 1;
  } else if (type == "kdtree2") {
    method = 2;
  } else if (type == "notree") {
    method = 0;
    bucketSize = dim(x)[2L];
  } else {
    stop("'type' must be 'kdtree0', 'kdtree1', 'kdtree2', or 'notree'");
  }

  if (dim(Xbal)[1L] != dim(Xspread)[2L])
    stop("the size of 'Xbal' and 'Xspread' does not match");

  if (length(prob == 1))
    stop("'prob' must be a vector of probabilities");

  if (length(prob) != dim(Xbal)[1L])
    stop("the size of 'prob' and 'Xbal' does not match");

  if (eps < 0.0 || 1e-4 < eps)
    stop("'eps' must be in [0.0, 1e-4]");

  result = .lcube_cpp(prob, Xbal, Xspread, bucketSize, method, eps);

  return(result);
}
