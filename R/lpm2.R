# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Local Pivotal Method 2
#'
#' @inherit lpm1 description details sections params return references
#'
#' @param prob A vector of length N with inclusion probabilities, or an integer > 1.
#' If an integer n, then the sample will be drawn with equal probabilities n / N.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = lpm2(prob, x);
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
#'   s = lpm2(prob, x);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
lpm2 = function(
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
    if (prob %% 1 != 0 || prob < 1)
      stop("'prob' must be a vector of probabilities or a single integer > 0");

    result = .lpm2_int_cpp(prob, x, bucketSize, method);
  } else {
    if (length(prob) != dim(x)[2L])
      stop("the size of 'prob' and 'x' does not match");

    if (eps < 0.0 || 1e-4 < eps)
      stop("'eps' must be in [0.0, 1e-4]");

    result = .lpm2_cpp(prob, x, bucketSize, method, eps);
  }

  return(result);
}
