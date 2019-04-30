#' Calculates geometric mean
#'
#' @param x - input vector
#' @param na.rm - should NA values be removed?
#'
#' @return
#' @export
#'
#' @examples
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

