# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Spatial balance
#'
#' @family measure
#'
#' @description
#' Calculates the spatial balance of a sample.
#'
#' @details
#' About voronoi and sumofsquares
#'
#' @templateVar xspread x
#' @templateVar integerprob TRUE
#' @template kdtrees_template
#' @template x_template
#' @template probs_template
#'
#' @param sample A vector of sample indices.
#' @param measure The type of balance measure to use.
#' Must be one of `"voronoi"`, `"sumofsquares"`.
#'
#' @return The balance measure of the provided sample.
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

  N = dim(x)[2L];
  method = .kdtree_method_check(type, bucketSize);
  bucketSize = .kdtree_bucket_check(N, type, bucketSize);
  prob = .prob_expand(prob, N);

  if (N < length(sample))
    stop("'sample' must be a vector of unique indices");

  if (measure == "voronoi") {
    result = .sb_voronoi_cpp(prob, x, sample, bucketSize, method);
  } else if (measure == "sumofsquares") {
    result = .sb_sumofsquares_cpp(x, sample, bucketSize, method);
  } else {
    stop("'balance' must be 'voronoi', or 'sumofsquares'");
  }

  return(result);
}
