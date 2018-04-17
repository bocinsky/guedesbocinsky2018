#' Modulate a record by a deviation of its standard score.
#'
#' @param data.df A dataframe of temperature data
#' @param temp.change.sd The change in standard score (z-score) of the data
#' @param t.base A base temperature for calculating GDD
#' @param t.cap An optional cap temperature for calculating GDD
#'
#' @return The modulated GDD
#' @export
sdModulator <- function(data.df, temp.change.sd, t.base, t.cap = NULL) {
  tmin <- ((data.df[["TMIN"]][, "SD"] / 10) * temp.change.sd) + (data.df[["TMIN"]][, "MEAN"] / 10)
  tmax <- ((data.df[["TMAX"]][, "SD"] / 10) * temp.change.sd) + (data.df[["TMAX"]][, "MEAN"] / 10)

  return(sum(calcGDD(tmin, tmax, t.base = t.base, t.cap = t.cap)))
}
