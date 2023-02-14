# **********************************************
# Author: Wilmer Prentius
# Licence: GPL (>=2)
# **********************************************

lpm1 = function(
  prob,
  x,
  type = "kdtree2",
  bucketSize = 40
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

    result = .lpm1_int_cpp(prob, x, bucketSize, method);
  } else {
    if (length(prob) != dim(x)[2L])
      stop("the size of 'prob' and 'x' does not match");

    result = .lpm1_cpp(prob, x, bucketSize, method);
  }

  return(result);
}
