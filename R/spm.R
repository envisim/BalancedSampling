# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

#' Sequential Pivotal Method
#'
#' @description
#' Selects samples with prescribed inclusion probabilities from a finite population.
#'
#' @details
#' If \code{prob} sum to an integer n, a fixed sized sample (n) will be produced.
#' The resulting samples are well spread in the list, similar to systematic sampling.
#' In each of the, at most, N steps, the two undecided units with the smallest
#' indices are selected to compete against each other.
#'
#' @template sampling_template
#' @template probs_template
#'
#' @references
#' Deville, J.-C. and  Tillé, Y. (1998).
#' Unequal probability sampling without replacement through a splitting method.
#' Biometrika 85, 89-101.
#'
#' Chauvet, G. (2012).
#' On a characterization of ordered pivotal sampling.
#' Bernoulli, 18(4), 1320-1340.
#'
#' @examples
#' \dontrun{
#' set.seed(12345);
#' N = 1000;
#' n = 100;
#' prob = rep(n/N, N);
#' x = matrix(runif(N * 2), ncol = 2);
#' s = spm(prob);
#' plot(x[, 1], x[, 2]);
#' points(x[s, 1], x[s, 2], pch = 19);
#'
#' set.seed(12345);
#' prob = c(0.2, 0.25, 0.35, 0.4, 0.5, 0.5, 0.55, 0.65, 0.7, 0.9);
#' N = length(prob);
#' ep = rep(0L, N);
#' r = 10000L;
#' for (i in seq_len(r)) {
#'   s = spm(prob);
#'   ep[s] = ep[s] + 1L;
#' }
#' print(ep / r);
#' }
#'
spm = function(
  prob,
  eps = 1e-12
) {
  if (length(prob) == 1)
    stop("'prob' must be a vector of probabilities");

  prob = as.numeric(prob);
  .eps_check(eps);

  result = .spm_cpp(prob, eps);

  return(result);
}
