# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Hierarchical Local Pivotal Method 2
#'
#' @inherit lpm1 details sections params return references
#'
#' @description
#' Selects an initial sample using the [lpm2()], and then splits this sample into
#' subsamples of given \code{sizes} using successive, hierarchical selection with
#' the [lpm2()].
#' The method is used to select several subsamples, such that each subsample, and
#' the combination (i.e. the union of all subsamples), is spatially balanced.
#'
#' @section Inclusion probabilities:
#' The inclusion probabilities \code{prob} _must_ sum to an integer n.
#' The sizes of the subsamples \code{sum(sizes)} _must_ sum to the same integer n.
#'
#' @param sizes A vector of integers containing the sizes of the subsamples.
#' \code{sum(sizes) = sum(prob)} must hold.
#'
#' @return A matrix with the population indices of the combined sample in the
#' first column, and the associated subsample in the second column.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' sizes = c(10, 20, 30, 40);
#' s = hlpm2(prob, x, sizes);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#' }
#'
hlpm2 = function(
  prob,
  x,
  sizes,
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
    bucketSize = dim(x)[2L];
  } else {
    stop("'type' must be 'kdtree0', 'kdtree1', 'kdtree2', or 'notree'");
  }

  if (bucketSize %% 1 != 0 || bucketSize < 1)
    stop("'bucketSize' must be integer > 0");

  if (eps < 0.0 || 1e-4 < eps)
    stop("'eps' must be in [0.0, 1e-4]");

  if (length(prob) == 1) {
    if (prob %% 1 != 0 || prob < 1 || prob > dim(x)[2L])
      stop("'prob' must be a vector of probabilities or a single integer in [0, N]");

    prob = rep(prob / dim(x)[2L], dim(x)[2L]);
  } else {
    if (length(prob) != dim(x)[2L])
      stop("the size of 'prob' and 'x' does not match");
  }

  probsum = sum(prob);

  if (probsum %% 1 != 0)
    stop("'prob' must sum to an integer");

  if (probsum != sum(sizes))
    stop("'sizes' must sum to an integer same as the sum of 'prob'");

  result = .hlpm2_cpp(prob, x, sizes, bucketSize, method, eps);

  return(result);
}
