# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Stratified doubly balanced sampling with pooling of landing phases
#'
#' @inherit lpm1 params return
#' @inherit cubestratified params
#' @inherit lcube params details
#' @inheritSection cube Inclusion probabilities
#' @inheritSection lpm1 k-d-trees
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
#' strata = c(rep(1L, 100), rep(2L, 200), rep(3L, 300), rep(4L, 400));
#' s = lcubestratified(prob, x, prob, strata);
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
#'   s = lcubestratified(prob, x, prob, strata);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
lcubestratified = function(
  prob,
  Xspread,
  Xbal,
  strata,
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

  if (length(strata) != dim(Xbal)[1L])
    stop("the size of 'strata' and 'Xbal' does not match");

  if (dim(Xbal)[1L] != dim(Xspread)[2L])
    stop("the size of 'Xbal' and 'Xspread' does not match");

  if (length(prob == 1))
    stop("'prob' must be a vector of probabilities");

  if (length(prob) != dim(Xbal)[1L])
    stop("the size of 'prob' and 'Xbal' does not match");

  if (eps < 0.0 || 1e-4 < eps)
    stop("'eps' must be in [0.0, 1e-4]");

  result = .lcube_stratified_cpp(prob, Xbal, Xspread, strata, bucketSize, method, eps);

  return(result);
}
