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

  assertthat::assert_that(is.vector(x),msg = "Please provide numeric vector as an input for gm_mean")
  assertthat::assert_that(length(x)>0,msg = "Empty vector provided as input")
  assertthat::assert_that(assertive::is_a_bool(na.rm),msg = "Please provide boolean value for na.rm option")
  assertthat::assert_that(assertive::is_numeric(x),msg = "Please provide numeric vector as input")
  if (length(x)==1) {
    gm_mean=x[1]
  }
  else {
    gm_mean = exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  return(gm_mean)
}

