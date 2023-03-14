# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Local Pivotal Method 1
#'
#' @description
#' Selects spatially balanced samples with prescribed inclusion probabilities
#' from a finite population using the Local Pivotal Method 1 (LPM1).
#'
#' @details
#' If \code{prob} sum to an integer n, a fixed sized sample (n) will be produced.
#'
#' @templateVar xspread x
#' @templateVar integerprob TRUE
#' @template sampling_template
#' @template kdtrees_template
#' @template x_template
#' @template probs_template
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

  N = dim(x)[2L];
  method = .kdtree_method_check(type, bucketSize);
  bucketSize = .kdtree_bucket_check(N, type, bucketSize);
  .eps_check(eps);

  if (.prob_integer_test(prob, N)) {
    result = .lpm_int_cpp(1, prob, x, bucketSize, method);
  } else {
    prob = .prob_check(prob, N);
    result = .lpm_cpp(1, prob, x, bucketSize, method, eps);
  }

  return(result);
}
