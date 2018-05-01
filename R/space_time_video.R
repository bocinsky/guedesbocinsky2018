#' A space-time video
#'
#' Given a raster brick (and possibly uncertainties/errors),
#' plot a video of the brick through time, and a
#' graph of the average across space, with a bar indicating year.
#'
#' @param the_brick The data to plot
#' @param the_brick_upper The upper CI data
#' @param the_brick_lower The lower CI data
#' @param out_file Where to write the output
#' @param title The Title of the plot
#' @param time The timeseries of the plot
#' @param timelim The time limits of the plot
#' @param timeaxis The time limits of the axis
#' @param timelab The label of the time axis
#' @param zbreaks The color breaks for plotting
#' @param zlab The label for the data
#' @param zaxis The axis breaks for the data
#' @param zcolors The colors for the data
#' @param fig_width The width of the figure
#' @param graph_height The height of the time series graph
#' @param margin The margin
#' @param pt_size The point size
#' @param smooth Should the timeseries be smoothed by kernel estimation
#' @param length The length (in seconds) of the video
#' @param extra_plot_fun Any extra plotting functions
#' @param extra_legend_fun Any extra legend functions
#'
#' @return NULL
#' @keywords internal
space_time_video <- function(the_brick,
                             the_brick_upper = NULL,
                             the_brick_lower = NULL,
                             out_file,
                             line = 75,
                             title = NULL,
                             time,
                             timelim,
                             timeaxis,
                             timelab = "Year AD",
                             zbreaks = NULL,
                             zlab,
                             zaxis,
                             zcolors = rev(RColorBrewer::brewer.pal(11, "RdYlBu")),
                             fig_width = 6.5,
                             graph_height = 1.5,
                             margin = 0.1,
                             pt_size = 8,
                             smooth = FALSE,
                             length = 60, # Length in seconds of video
                             extra_plot_fun = NULL,
                             extra_legend_fun = NULL) {
  out_dir <- tempfile()
  dir.create(out_dir)

  if (isTRUE(smooth)) {
    smoother <- stats::dnorm(seq(-10, 10, 1), sd = 5)
  } else if (is.numeric(smooth)) {
    smoother <- smooth
  } else {
    smoother <- NULL
  }

  mean.all <- mean(the_brick[], na.rm = T)
  mean.spatial <- raster::mean(the_brick, na.rm = T)
  mean.temporal <- raster::cellStats(the_brick, mean, na.rm = T)
  ci.temporal <- raster::cellStats(the_brick, stats::quantile, probs = c(0.25, 0.75), na.rm = T)

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
    mean.spatial.lower <- raster::mean(the_brick_lower, na.rm = T)
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

  for (layer in 1:raster::nlayers(the_brick)) {
    grDevices::png(
      filename = stringr::str_c(out_dir, "/image",
                                stringi::stri_pad_left(layer, width = 4, pad = 0), ".png"),
      width = fig_width,
      height = fig_height,
      units = "in",
      type = "cairo-png",
      antialias = "none",
      pointsize = 8,
      res = 600
    )

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

    raster::plot(the_brick[[layer]],
      maxpixels = raster::ncell(the_brick),
      zlim = range(zbreaks),
      breaks = zbreaks,
      add = T,
      col = colors,
      colNA = "gray90",
      useRaster = TRUE,
      legend = FALSE
    )

    if ((min(the_brick[[layer]][], na.rm = TRUE) != max(the_brick[[layer]][], na.rm = TRUE))) {
      raster::contour(the_brick[[layer]],
        maxpixels = raster::ncell(the_brick),
        levels = line,
        drawlabels = FALSE,
        col = "white",
        lwd = 1.25,
        add = T
      )
    }

    if (!is.null(the_brick_upper) & (min(the_brick_upper[[layer]][], na.rm = TRUE) != max(the_brick_upper[[layer]][], na.rm = TRUE))) {
      raster::contour(the_brick_upper[[layer]],
        maxpixels = raster::ncell(the_brick_upper),
        levels = line,
        drawlabels = FALSE,
        col = "white",
        lwd = 0.75,
        lty = 1,
        add = T
      )
    }

    if (!is.null(the_brick_lower) & (min(the_brick_lower[[layer]][], na.rm = TRUE) != max(the_brick_lower[[layer]][], na.rm = TRUE))) {
      raster::contour(the_brick_lower[[layer]],
        maxpixels = raster::ncell(the_brick_lower),
        levels = line,
        drawlabels = FALSE,
        col = "white",
        lwd = 0.75,
        lty = 1,
        add = T
      )
    }

    if (!is.null(extra_plot_fun)) {
      extra_plot_fun(years = time[[layer]])
    }

    if (!is.null(extra_legend_fun)) {
      extra_legend_fun()
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
      h = line,
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

    graphics::text(
      x = fig_width,
      y = max(zbreaks),
      labels = time[[layer]] %>%
        stringr::str_c(" Years BP"),
      adj = c(1, -1),
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

    graphics::abline(
      v = time[[layer]],
      lty = 1,
      lwd = 1.5,
      col = "#cb181d",
      xpd = FALSE
    )

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
  }

  movie.width <- 1600
  movie.height <- round(movie.width * (fig_height / fig_width))
  if (movie.height %% 2 != 0) movie.height <- movie.height + 1

  # Create the video
  fps <- raster::nlayers(the_brick) / length
  system(stringr::str_c("ffmpeg",
                        " -r ", fps,
                        " -i ", out_dir, "/image%04d.png",
                        " -s:v ", movie.width, "x", movie.height,
                        " -c:v libx264",
                        " -crf 30",
                        " -profile:v High",
                        " -pix_fmt yuv420p",
                        " ",
                        out_file, " -y"))
}
