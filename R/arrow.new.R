utils::globalVariables(c("inch.x",
                         "inch.y",
                         "plot.extent"))

#' A north arrow.
#'
#' @param x The x position of the arrow
#' @param y The y position of the arrow
#' @param height The height of the arrow
#' @param width The width of the arrow
#' @param adj The adjustment for the arrow
#'
#' @export
arrow.new <- function(x, y, height, width, adj) {
  graphics::arrows(
    x0 = raster::xmin(plot.extent) + 0.2 * inch.x,
    y0 = raster::ymin(plot.extent) + 0.3 * inch.y,
    x1 = raster::xmin(plot.extent) + 0.2 * inch.x,
    y1 = raster::ymin(plot.extent) + 0.6 * inch.y,
    length = 0.1,
    lwd = 1.5,
    lend = 1
  )
  graphics::text(
    labels = "N",
    x = raster::xmin(plot.extent) + 0.2 * inch.x,
    y = raster::ymin(plot.extent) + 0.3 * inch.y,
    adj = c(0.5, 0),
    cex = 1.5,
    font = 2
  )
}
