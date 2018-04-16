climate_plotter <- function(data, station, element){
  weather <- data$weather[[station]][[element]] %>%
    as_data_frame()
  
  weather <- data.frame(DOY=as.numeric(strftime(as.POSIXlt(paste(rep(weather$YEAR,each=31),rep(weather$MONTH,each=31),rep(1:31,times=nrow(weather)), sep='.'), format="%Y.%m.%d"), format = "%j")), DATA=as.numeric(t(weather[,-1:-2])) )
  weather <- weather[!is.na(weather$DATA) & !is.na(weather$DOY),]
  
  ymin <- min(weather$DATA/10) %>% round_any(10, f = floor)
  ymax <- max(weather$DATA/10) %>% round_any(10, f = ceiling)
  
  climate <- data$climatology[[station]][[element]] %>%
    as_data_frame()
  
  p <- ggplot() +
    geom_point(data = weather, mapping = aes(x = DOY, y = DATA/10), size = 0.1) +
    geom_line(data = climate, mapping = aes(x = DOY, y = MEAN/10), color = "red") +
    geom_line(data = climate, mapping = aes(x = DOY, y = (MEAN + SD)/10), color = "dodgerblue") +
    geom_line(data = climate, mapping = aes(x = DOY, y = (MEAN - SD)/10), color = "dodgerblue") +
    coord_polar(theta = "x") +
    ylim(ymin,ymax) +
    xlim(1,366) +
    scale_x_continuous(breaks=c(1,32,60,92,122,153,184,214,245,275,306,336),
                       labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")) +
    theme_linedraw() +
    theme(axis.title.x=element_blank())
  
  if(element == "TMIN"){
    p <- p +
      ylab("Minimum temperature (ºC)")
  }else if(element == "TMAX"){
    p <- p +
      ylab("Maximum temperature (ºC)")
  }else if(element == "PRCP"){
    p <- p +
      ylab("Precipitation (mm)")
  }
  p
}