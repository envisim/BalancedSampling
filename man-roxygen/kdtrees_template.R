#' @section k-d-trees:
#' The \code{type}s "kdtree" creates k-d-trees with terminal node bucket sizes
#' according to \code{bucketSize}.
#'
#' \itemize{
#'   \item{"kdtree0"} creates a k-d-tree using a median split on alternating variables.
#'   \item{"kdtree1"} creates a k-d-tree using a median split on the largest range.
#'   \item{"kdtree2"} creates a k-d-tree using a sliding-midpoint split.
#'   \item{"notree"} does a naive search for the nearest neighbour.
#' }
#'
#' @param type The method used in finding nearest neighbours.
#' Must be one of \code{"kdtree0"}, \code{"kdtree1"}, \code{"kdtree2"}, and
#' \code{"notree"}.
#' @param bucketSize The maximum size of the terminal nodes in the k-d-trees.
#'
