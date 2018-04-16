# This function smooths periodic (annual) data in the GHCN format.
# First, it converts all data to DOY data (366-day-year).
# Then, it calculates a periodic generalized additive model 
# from the data using cyclic cubic regression splines for smoothing.
# It then calculates the daily-wise standard deviations of the residuals,
# and again smooths using the same type of GAM.
# Optionally, you may plot the results.
calcDailyMeanSD <- function(daily.data, unwrapped=F, plot=F){
  if(!unwrapped){
    # Unwrap data and convert to POSIX Day of Year (DoY)
    daily.data.unwrap <- data.frame(DOY=as.numeric(strftime(as.POSIXlt(paste(rep(daily.data$YEAR,each=31),rep(daily.data$MONTH,each=31),rep(1:31,times=nrow(daily.data)), sep='.'), format="%Y.%m.%d"), format = "%j")), DATA=as.numeric(t(daily.data[,-1:-2])) )
    daily.data.unwrap <- daily.data.unwrap[!is.na(daily.data.unwrap$DATA) & !is.na(daily.data.unwrap$DOY),]
  }else{
    daily.data.unwrap <- daily.data
  }

  # fit a generalized additive model to data as smooth, periodic function of DoY
  mean.model <- gam(DATA ~ s(DOY, bs='cc'), data=daily.data.unwrap)
  
  # generate predictions from our DoY model
  mean.predictions <- predict(mean.model, data.frame(DOY = 1:366))

  # Get the standard deviation of the residuals for each day
  daily.data.unwrap$RES <- residuals(mean.model)
  residual.sds <- aggregate(daily.data.unwrap$RES,by=list(daily.data.unwrap$DOY),sd)
  names(residual.sds) <- c("DOY","SD")
  sd.model <- gam(SD ~ s(DOY, bs='cc'), data=residual.sds)
  sd.predictions <- predict(sd.model, data.frame(DOY = 1:366))
  
  if(plot){
    polar.plot(lengths=daily.data.unwrap$DATA, polar.pos=(daily.data.unwrap$DOY)*360/365, labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov", "Dec"), label.pos=c(0,32,60,92,122,153,184,214,245,275,306,336)*360/365, rp.type='s', clockwise=TRUE, start=1, radial.lim=range(daily.data.unwrap$DATA), cex=0.1)
    polar.plot(lengths=as.numeric(mean.predictions), polar.pos=as.numeric(names(mean.predictions))*360/366, rp.type='p', clockwise=TRUE, start=0, radial.lim=range(daily.data.unwrap$DATA), lwd=3, line.col='dodgerblue', add=TRUE)
    polar.plot(lengths=as.numeric(mean.predictions)+as.numeric(sd.predictions), polar.pos=as.numeric(names(mean.predictions))*360/366, rp.type='p', clockwise=TRUE, start=0, add=TRUE, lty=1, radial.lim=range(daily.data.unwrap$DATA), line.col='red', lwd=2)
    polar.plot(lengths=as.numeric(mean.predictions)-as.numeric(sd.predictions), polar.pos=as.numeric(names(mean.predictions))*360/366, rp.type='p', clockwise=TRUE, start=0, add=TRUE, lty=1, radial.lim=range(daily.data.unwrap$DATA), line.col='red', lwd=2)
#     polar.plot(lengths=rep(mean(as.numeric(mean.predictions)),366), polar.pos=as.numeric(names(mean.predictions))*360/366, rp.type='p', clockwise=TRUE, start=0, radial.lim=range(daily.data.unwrap$DATA), lwd=2, line.col='dodgerblue', add=TRUE)
  }
  
  return(do.call(cbind,list(DOY=1:366, MEAN=mean.predictions, SD=sd.predictions)))
}