#' A function to calculate GDD, given daily maximum and minimum temperature series,
#' and maximum and minimum growth temperatures.
#'
#' GDD is calculated as
#' $[(TMINc+TMAXc)/2]-TBASE$
#' Where TMINc and TMAXc are measured temperatures bounded by TBASE as a minimum
#' and perhaps a TCAP as a maximum.
#'
#' @param tmin.vector
#' @param tmax.vector
#' @param t.base
#' @param t.cap
#'
#' @return
#' @export
#'
#' @examples
calcGDD <- function(tmin.vector,
                    tmax.vector,
                    t.base,
                    t.cap = NULL) {
  # Floor all temps at t.base
  tmin.vector[tmin.vector < t.base] <- t.base
  tmax.vector[tmax.vector < t.base] <- t.base

  # Ceiling all temps at t.cap
  if (!is.null(t.cap)) {
    tmin.vector[tmin.vector > t.cap] <- t.cap
    tmax.vector[tmax.vector > t.cap] <- t.cap
  }

  # GDD is t.avg-t.base
  t.avg <- (tmin.vector + tmax.vector) / 2

  GDD <- t.avg - t.base

  return(GDD)
}
