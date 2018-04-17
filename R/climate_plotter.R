utils::globalVariables(c("DOY",
                         "DATA",
                         "MEAN",
                         "SD"))
#' Plot climatology for a GHCN station and element (variable).
#'
#' @param data GHCN data
#' @param station A station name
#' @param element An element (variable) name
#'
#' @return A ggplot2 object
#' @importFrom magrittr %>%
#' @export
climate_plotter <- function(data, station, element) {
  weather <- data$weather[[station]][[element]] %>%
    tibble::as_data_frame()

  weather <- data.frame(DOY = as.numeric(strftime(as.POSIXlt(paste(rep(weather$YEAR, each = 31),
                                                                   rep(weather$MONTH, each = 31),
                                                                   rep(1:31, times = nrow(weather)),
                                                                   sep = "."),
                                                             format = "%Y.%m.%d"),
                                                  format = "%j")),
                        DATA = as.numeric(t(weather[, -1:-2])))
  weather <- weather[!is.na(weather$DATA) & !is.na(weather$DOY), ]

  ymin <- min(weather$DATA / 10) %>% plyr::round_any(10, f = floor)
  ymax <- max(weather$DATA / 10) %>% plyr::round_any(10, f = ceiling)

  climate <- data$climatology[[station]][[element]] %>%
    tibble::as_data_frame()

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = weather, mapping = ggplot2::aes(x = DOY, y = DATA / 10), size = 0.1) +
    ggplot2::geom_line(data = climate, mapping = ggplot2::aes(x = DOY, y = MEAN / 10), color = "red") +
    ggplot2::geom_line(data = climate, mapping = ggplot2::aes(x = DOY, y = (MEAN + SD) / 10), color = "dodgerblue") +
    ggplot2::geom_line(data = climate, mapping = ggplot2::aes(x = DOY, y = (MEAN - SD) / 10), color = "dodgerblue") +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::ylim(ymin, ymax) +
    ggplot2::xlim(1, 366) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 32, 60, 92, 122, 153, 184, 214, 245, 275, 306, 336),
      labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    ) +
    ggplot2::theme_linedraw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  if (element == "TMIN") {
    p <- p +
      ggplot2::ylab("Minimum temperature (\u00B0C)")
  } else if (element == "TMAX") {
    p <- p +
      ggplot2::ylab("Maximum temperature (\u00B0C)")
  } else if (element == "PRCP") {
    p <- p +
      ggplot2::ylab("Precipitation (mm)")
  }
  p
}
