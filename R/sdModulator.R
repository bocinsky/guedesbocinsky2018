#' Modulate a record by a deviation of its standard score.
#'
#' @param data.df
#' @param temp.change.sd
#' @param t.base
#' @param t.cap
#'
#' @return
#' @export
#'
#' @examples
sdModulator <- function(data.df, temp.change.sd, t.base, t.cap = NULL){
  tmin <- ((data.df[['TMIN']][,"SD"]/10)*temp.change.sd) + (data.df[['TMIN']][,"MEAN"]/10)
  tmax <- ((data.df[['TMAX']][,"SD"]/10)*temp.change.sd) + (data.df[['TMAX']][,"MEAN"]/10)

  return(sum(calcGDD(tmin, tmax, t.base=t.base, t.cap=t.cap)))
}
