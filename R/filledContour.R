#' Add a filled contour to a plot.
#'
#' @param x a raster object
#' @param y Integer. The layer number of x (if x has multiple layers)
#' @param maxpixels the maximum number of pixels to plot
#' @param levels The levels to plot
#' @param col a vector of colors
#'
#' @return An additional plot on the current graphics device
#' @export
filledContourAdd <- function(x,
                             y = 1,
                             maxpixels = 1e+05,
                             levels,
                             col) {
  if (raster::nlayers(x) > 1) {
    y <- min(max(1, y), raster::nlayers(x))
    x <- raster::raster(x, y)
  }
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE, useGDAL = TRUE)
  X <- raster::xFromCol(x, 1:ncol(x))
  Y <- raster::yFromRow(x, nrow(x):1)
  Z <- t(matrix(raster::getValues(x), ncol = x@ncols, byrow = TRUE)[raster::nrow(x):1, ])
  graphics::.filled.contour(x = X, y = Y, z = Z, levels = levels, col = col)
}
