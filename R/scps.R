# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Spatially Correlated Poisson Sampling
#'
#' @description
#' Selects spatially balanced samples with prescribed inclusion probabilities
#' from a finite population using Spatially Correlated Poisson Sampling (SCPS).
#'
#' @details
#' If `prob` sum to an integer n, a fixed sized sample (n) will be produced.
#' The implementation uses the maximal weight strategy, as specified in
#' Grafström (2012).
#'
#' # Coordinated SCPS
#' If `rand` is supplied, coordinated SCPS will be performed.
#' The algorithm for coordinated SCPS differs from the SCPS algorithm, as
#' uncoordinated SCPS chooses a unit to update randomly, whereas coordinated SCPS
#' traverses the units in the supplied order.
#' This has a small impact on the efficiency of the algorithm for coordinated SCPS.
#'
#' @templateVar xspread x
#' @templateVar integerprob TRUE
#' @template sampling_template
#' @template kdtrees_template
#' @template x_template
#' @template probs_template
#'
#' @param rand A vector of length N with random numbers.
#' If this is supplied, the decision of each unit is taken with the corresponding
#' random number. This makes it possible to coordinate the samples.
#'
#' @references
#' Grafström, A. (2012).
#' Spatially correlated Poisson sampling.
#' Journal of Statistical Planning and Inference, 142(1), 139-147.
#'
#' Friedman, J. H., Bentley, J. L., & Finkel, R. A. (1977).
#' An algorithm for finding best matches in logarithmic expected time.
#' ACM Transactions on Mathematical Software (TOMS), 3(3), 209-226.
#'
#' Maneewongvatana, S., & Mount, D. M. (1999, December).
#' It’s okay to be skinny, if your friends are fat.
#' In Center for geometric computing 4th annual workshop on computational geometry (Vol. 2, pp. 1-8).
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = scps(prob, x);
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
#'   s = scps(prob, x);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
scps = function(
  prob,
  x,
  rand = NULL,
  type = "kdtree2",
  bucketSize = 50,
  eps = 1e-12
) {
  if (!is.matrix(x)) {
    x = t(as.matrix(x));
  } else {
    x = t(x);
  }

  N = dim(x)[2L];
  method = .kdtree_method_check(type, bucketSize);
  bucketSize = .kdtree_bucket_check(N, type, bucketSize);
  .eps_check(eps);
  prob = .prob_expand(prob, N);

  if (rand != FALSE && is.vector(rand)) {
    if (length(rand) != N)
      stop("the size of 'rand' and 'x' does not match");

    result = .scps_coord_cpp(prob, x, rand, bucketSize, method, eps);
  } else {
    result = .scps_cpp(prob, x, bucketSize, method, eps);
  }

  return(result);
}
