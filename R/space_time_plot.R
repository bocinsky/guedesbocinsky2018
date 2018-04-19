#' A space-time plot.
#'
#' Given a raster brick (and possibly uncertainties/errors),
#' plot a map of the average across time, and a
#' graph of the average across space.
#'
#' @param the_brick
#' @param the_brick_upper
#' @param the_brick_lower
#' @param out_file
#' @param title
#' @param time
#' @param timelim
#' @param timeaxis
#' @param timelab
#' @param zbreaks
#' @param zlab
#' @param zaxis
#' @param zcolors
#' @param fig_width
#' @param graph_height
#' @param margin
#' @param pt_size
#' @param smooth
#' @param extra_plot_fun
#' @param extra_legend_fun
#'
#' @return
#' @keywords internal
#' @importFrom magrittr %<>%
space_time_plot <- function(the_brick,
                            the_brick_upper = NULL,
                            the_brick_lower = NULL,
                            out_file,
                            title = NULL,
                            time,
                            timelim,
                            timeaxis,
                            timelab = "Year AD",
                            zbreaks = NULL,
                            # zbreaks_mid_range = NULL,
                            zlab,
                            zaxis,
                            zcolors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                            # zlim_colors_skew = FALSE,
                            fig_width = 6.5,
                            graph_height = 1.5,
                            margin = 0.1,
                            pt_size = 8,
                            smooth = FALSE,
                            extra_plot_fun = NULL,
                            extra_legend_fun = NULL) {
  if (isTRUE(smooth)) {
    smoother <- stats::dnorm(seq(-10, 10, 1), sd = 5)
  } else if (is.numeric(smooth)) {
    smoother <- smooth
  } else {
    smoother <- NULL
  }

  if (!raster::inMemory(the_brick)) {
    the_brick %<>%
      raster::readAll()
  }

  mean.all <- mean(the_brick[], na.rm = T)
  mean.spatial <- raster::mean(the_brick, na.rm = T)
  mean.temporal <- raster::cellStats(the_brick, mean, na.rm = T)
  ci.temporal <- raster::cellStats(the_brick,
                                   stats::quantile,
                                   probs = c(0.25, 0.75),
                                   na.rm = T)

  if (!is.null(smoother)) {
    mean.temporal %<>% stats::filter(filter = smoother)
  }

  if (!is.null(the_brick_upper)) {
    mean.all.upper <- mean(the_brick_upper[], na.rm = T)
    mean.spatial.upper <- raster::mean(the_brick_upper, na.rm = T)
    mean.temporal.upper <- raster::cellStats(the_brick_upper, mean, na.rm = T)
  }

  if (!is.null(the_brick_lower)) {
    mean.all.lower <- mean(the_brick_lower[], na.rm = T)
    mean.spatial.lower <-  raster::mean(the_brick_lower, na.rm = T)
    mean.temporal.lower <- raster::cellStats(the_brick_lower, mean, na.rm = T)
  }

  ym <- mean(c(the_brick@extent@ymax, the_brick@extent@ymin))

  aspect <- ifelse(raster::isLonLat(the_brick),
    (ncol(the_brick) / nrow(the_brick)) / (1 / cos(ym * pi / 180)),
    nrow(the_brick) / ncol(the_brick)
  )

  plot_width <- fig_width - (margin * 2)
  plot_height <- plot_width / aspect
  fig_height <- plot_height + (margin * 3) + graph_height

  if (is.null(zbreaks)) {
    zbreaks <- seq(min(the_brick[], na.rm = T),
      max(the_brick[], na.rm = T),
      length.out = 100
    )
  }

  colors <- grDevices::colorRampPalette(zcolors)(length(zbreaks) - 1)

  grDevices::cairo_pdf(
    filename = out_file,
    width = fig_width,
    height = fig_height,
    antialias = "none",
    pointsize = 8,
    fallback_resolution = 600
  )

  # grDevices::pdf(
  #   file = out_file,
  #   width = fig_width,
  #   height = fig_height,
  #   pointsize = 8
  # )

  graphics::par(
    mai = c(
      graph_height + (margin * 2),
      margin,
      margin,
      margin
    ),
    xpd = F
  )

  graphics::plot(1,
    type = "n",
    xlab = "",
    ylab = "",
    xlim = c(raster::extent(the_brick)@xmin, raster::extent(the_brick)@xmax),
    ylim = c(raster::extent(the_brick)@ymin, raster::extent(the_brick)@ymax),
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    main = ""
  )

  raster::plot(mean.spatial,
    maxpixels = raster::ncell(mean.spatial),
    zlim = range(zbreaks),
    breaks = zbreaks,
    add = T,
    col = colors,
    colNA = "gray90",
    useRaster = TRUE,
    legend = FALSE
  )

  raster::contour(mean.spatial,
    maxpixels = raster::ncell(mean.spatial),
    levels = 0.75,
    drawlabels = FALSE,
    col = "white",
    lwd = 1.25,
    add = T
  )

  if (!is.null(the_brick_upper)) {
    raster::contour(mean.spatial.upper,
      maxpixels = raster::ncell(mean.spatial.upper),
      levels = 0.75,
      drawlabels = FALSE,
      col = "white",
      lwd = 0.75,
      lty = 1,
      add = T
    )
  }

  if (!is.null(the_brick_lower)) {
    raster::contour(mean.spatial.lower,
      maxpixels = raster::ncell(mean.spatial.lower),
      levels = 0.75,
      drawlabels = FALSE,
      col = "white",
      lwd = 0.75,
      lty = 1,
      add = T
    )
  }

  if (!is.null(extra_plot_fun)) {
    extra_plot_fun()
  }

  graphics::par(mai = c(
    (margin * 2),
    margin,
    (margin * 3) + plot_height,
    margin
  ), xpd = T, new = T)

  graphics::plot(1,
    type = "n",
    xlab = "",
    ylab = "",
    xlim = c(0, fig_width),
    ylim = range(zbreaks),
    xaxs = "i",
    yaxs = "i",
    axes = FALSE,
    main = ""
  )

  legend.breaks <- seq(from = utils::head(zbreaks, 1),
                       to = utils::tail(zbreaks, 1),
                       length.out = (length(zbreaks) + 1))
  graphics::rect(
    col = colors,
    border = NA,
    ybottom = zbreaks[1:(length(zbreaks) - 1)],
    ytop = zbreaks[2:length(zbreaks)],
    xleft = 0.15,
    xright = 0.35,
    xpd = T
  )
  graphics::abline(
    h = 0.75,
    col = "white"
  )
  graphics::text(
    x = 0,
    y = mean(zbreaks),
    labels = zlab,
    adj = c(0.5, 1),
    cex = 0.9,
    srt = 90,
    font = 2
  )
  graphics::text(
    x = 0.5,
    y = c(
      utils::head(zbreaks, 1),
      utils::tail(zbreaks, 1),
      zaxis
    ),
    labels = c(
      utils::head(zbreaks, 1),
      utils::tail(zbreaks, 1),
      zaxis
    ),
    adj = c(0.5, 0.5),
    cex = 0.8
  )

  graphics::text(
    x = 0,
    y = max(zbreaks),
    labels = title,
    adj = c(0, -1),
    cex = 1,
    font = 2
  )

  graphics::par(
    mai = c(
      (margin * 2),
      margin * 8,
      (margin * 3) + plot_height,
      margin * 2
    ),
    xpd = T,
    new = T
  )
  graphics::plot(1,
                 type = "n",
                 xlab = "",
                 ylab = "",
                 xlim = timelim,
                 ylim = range(zbreaks),
                 xaxs = "i",
                 yaxs = "i",
                 axes = FALSE,
                 main = "")


  graphics::polygon(
    x = c(time, rev(time)),
    y = c(ci.temporal[1, ], rev(ci.temporal[2, ])),
    border = NA,
    col = "gray90"
  )
  graphics::lines(
    y = mean.temporal,
    x = time,
    lwd = 1.5
  )

  graphics::abline(
    h = mean.all,
    lty = 2,
    lwd = 1.5,
    xpd = FALSE
  )

  if (!is.null(the_brick_lower)) {
    graphics::lines(
      y = mean.temporal.lower,
      x = time,
      lwd = 0.5,
      lty = 1
    )

    graphics::abline(
      h = mean.all.lower,
      lty = 2,
      lwd = 0.5,
      xpd = FALSE
    )
  }

  if (!is.null(the_brick_upper)) {
    graphics::lines(
      y = mean.temporal.upper,
      x = time,
      lwd = 0.5,
      lty = 1
    )

    graphics::abline(
      h = mean.all.upper,
      lty = 2,
      lwd = 0.5,
      xpd = FALSE
    )
  }

  graphics::axis(2,
    at = c(
      utils::head(zbreaks, 1),
      utils::tail(zbreaks, 1),
      zaxis
    ),
    labels = F
  )

  graphics::par(
    mai = c(
      margin,
      margin * 8,
      (margin * 2) + plot_height,
      margin * 2
    ),
    xpd = T,
    new = T
  )
  graphics::plot(1,
                 type = "n",
                 xlab = "",
                 ylab = "",
                 xlim = timelim,
                 ylim = c(0, graph_height),
                 xaxs = "i",
                 yaxs = "i",
                 axes = FALSE,
                 main = "")
  graphics::segments(
    x0 = timeaxis,
    x1 = timeaxis,
    y0 = (margin * 1),
    y1 = graph_height - (margin * 1),
    col = "gray50",
    lty = 3
  )
  graphics::text(
    x = timeaxis,
    y = 0,
    labels = timeaxis,
    adj = c(0.5, 0),
    cex = 0.8
  )
  graphics::text(
    x = timelim[1],
    y = 0,
    labels = timelab,
    adj = c(-0.1, 0),
    cex = 0.9,
    font = 2
  )

  grDevices::dev.off()
  # distill(out_file)
}
