#' A simple scale bar.
#'
#' @param d the length of the scalebar in meters
#' @param xy the position of the scale bar
#' @param height the height of the scale bar in inches
#' @param line.offset the offset of the position of the scale bar
#' @param side which side of the plot the scale bar should be on
#' @param lab.side which side of the scale bar the label should be on
#' @param col the color of the scale bar
#' @param lonlat is the plot in georgaphic coordinates?
#' @param label The label
#' @param adj An adjustment for the scale bar
#' @param lwd the line width of the scale bar
#' @param ... Additional parameters passed on to text
#'
#' @return A scalebar is added to the current graphics device.
#' @export
scalebar.new <- function(d,
                         xy = NULL,
                         height = NULL,
                         line.offset=c(0, 0),
                         side="right",
                         lab.side="top",
                         col="black",
                         lonlat = NULL,
                         label,
                         adj = c(0.5, -0.5),
                         lwd = 2, ...) {
  pr <- graphics::par()
  if (is.null(lonlat)) {
    if (pr$usr[1] > -181 & pr$usr[2] < 181 & pr$yaxp[1] >
      -200 & pr$yaxp[2] < 200) {
      lonlat <- TRUE
    }
    else {
      lonlat <- FALSE
    }
  }

  if (lonlat) {
    lat <- mean(pr$yaxp[1:2])
    if (missing(d)) {
      dx <- (pr$usr[2] - pr$usr[1]) / 10
      d <- raster::pointDistance(
        cbind(0, lat), cbind(dx, lat),
        TRUE
      )
      d <- signif(d / 1000, 2)
      label <- NULL
    }
    p <- cbind(0, lat)
    dd <- raster:::.destPoint(p, d * 1000)
    dd <- dd[1, 1]
  } else {
    if (missing(d)) {
      d <- round(10 * (pr$usr[2] - pr$usr[1]) / 10) / 10
      label <- NULL
    }
    dd <- d
  }

  if (is.null(xy)) {
    padding <- c(5, 5) / 100
    parrange <- c(pr$usr[2] - pr$usr[1], pr$usr[4] - pr$usr[3])
    xy <- c(pr$usr[1] + (padding[1] * parrange[1]), pr$usr[3] +
      (padding[2] * parrange[2]))
  }

  xy <- xy + line.offset

  if (side == "right") {
    xstart <- xy[1]
    xend <- xy[1] + dd
  } else {
    xstart <- xy[1] - dd
    xend <- xy[1]
  }

  if (is.null(height)) {
    height <- dd * 0.1
  }

  graphics::rect(xleft = xstart, ybottom = xy[2], xright = xend, ytop = xy[2] + height, col = col, border = NA, lend = 1)

  #   lines(matrix(c(xstart, xy[2], xend, xy[2]), byrow = T, nrow = 2), lend=1, lwd = lwd, ...)

  if (missing(label)) {
    label <- paste(d)
  }
  if (is.null(label)) {
    label <- paste(d)
  }
  if (missing(adj)) {
    adj <- c(0.5, -0.2 - lwd / 20)
  }

  if (lab.side == "top") {
    graphics::text(mean(c(xstart, xend)), xy[2],
      labels = label, adj = c(0.5, -0.5), col = col,
      ...
    )
  } else if (lab.side == "right") {
    graphics::text(xend, xy[2],
      labels = label, adj = c(-0.1, 0), col = col,
      ...
    )
  } else if (lab.side == "left") {
    graphics::text(xstart, xy[2],
      labels = label, adj = c(1.1, 0), col = col,
      ...
    )
  } else if (lab.side == "bottom") {
    graphics::text(mean(c(xstart, xend)), xy[2],
      labels = label, adj = c(0.5, 1.5), col = col,
      ...
    )
  }
}
