# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Spatial balance
#'
#' Calculates the spatial balance of a sample.
#'
#' @inherit lpm1 params references
#'
#' @param prob A vector of length N with inclusion probabilities, or an integer > 1.
#' If an integer n, then the probabilities will be assumed to be equal n / N.
#' @param sample A vector of sample indices.
#' @param measure The type of balance measure to use.
#' Must be one of \code{"voronoi"}, \code{"sumofsquares"}.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 500;
#' n = 70;
#' prob = rep(n / N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = lpm2(prob, x);
#' b = sb(prob, x, s);
#' }

sb = function(
  prob,
  x,
  sample,
  measure = "voronoi",
  type = "kdtree2",
  bucketSize = 10
) {
  if (!is.matrix(x)) {
    x = t(as.matrix(x));
  } else {
    x = t(x);
  }

  if (dim(x)[2L] < length(sample))
    stop("'sample' must be a vector of unique indices");

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

    prob = rep(prob / dim(x)[2L], dim(x)[2L]);
  } else {
    if (length(prob) != dim(x)[2L])
      stop("the size of 'prob' and 'x' does not match");
  }

  if (measure == "voronoi") {
    result = .sb_voronoi_cpp(prob, x, sample, bucketSize, method);
  } else if (measure == "sumofsquares") {
    result = .sb_sumofsquares_cpp(x, sample, bucketSize, method);
  } else {
    stop("'balance' must be 'voronoi', or 'sumofsquares'");
  }

  return(result);
}
