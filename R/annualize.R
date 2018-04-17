## A function that sums a GHCN record over a set of months to generate annual estimates.
## Really, you could calculate gdd over any period in the year that you would like
annualize <- function(station, months=1:12){
  # first, remove unwanted months
  station <- station[station$MONTH %in% months,]
  
  # then, split by year
  station.years <- split(station,station$YEAR)
  
  # then, sum over each year
  years.out <- lapply(station.years,function(year){
    sum(year[,3:33],na.rm=T)
  })
  
  return(as.data.frame(t(unlist(years.out))))
}