#' Calculates geometric mean
#'
#' @param x input vector
#' @param na.rm should NA values be removed?
#'
#' @return geometric mean of values provided as an input
#' @export
#'
#' @examples
#' a <- rnorm(100,33,5)
#' gm_mean(a)
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

