## A function that cleans ghcn data by four iterated processes:
## (1) crop to year.range
## (2) purge years with gaps in the data longer than the cutoff length
## (3) interpolate across gaps <= the cutoff length
## (4) fill by average those days at the beginning and end of years with missing data
## (5) remove stations with less than min.years years of complete data
ghcnCleaner <- function(data.list, year.range, na.cutoff, min.years=length(year.range)){
  ## (1) crop each series to year.range
  data.list <- lapply(data.list,function(X){
    X <- X[X$YEAR %in% year.range,]
    return(X)
  })
  
  ## (2) purge years with gaps in the data longer than the cutoff length
  data.list <- lapply(data.list,function(X){
    # Split records by year
    all.records.list <- split(X,X$YEAR)
    
    # For each year, do the following
    complete.records.list <- lapply(all.records.list,function(annual.data){
      # If less than 12 months, return null
      if(nrow(annual.data)<12){
        return(NULL)
      }
      
      # Get just the climate info
      annual.records <- as.matrix(annual.data[,4:34])
      
      # Get the number of days per month in the records
      n.days <- monthDays(as.Date(paste(annual.data$YEAR,annual.data$MONTH,'01',sep='-')))
      
      ## Unnwrap each row, accounting for number of days in the month
      annual.records.unwrapped <- unwrapRows(annual.records,n.days)
      
      ## count consecutive missing values (these are the same)
      annual.records.rle <- rle(is.na(annual.records.unwrapped))
      
      ## If any of the annual spans of NAs are > the cutoff, purge that year
      if(any(annual.records.rle$lengths[annual.records.rle$values] > na.cutoff)){
        return(NULL)
      }else{
        ## Otherwise, return the year
        return(annual.data)
      } 
    })
    
    complete.records.list <- complete.records.list[!sapply(complete.records.list, is.null)]
    
    if(length(complete.records.list) == 0 || length(complete.records.list)<min.years){
      return(NULL)
    }
    
    # rbind all complete years back together for interpolation
    complete.records.list <- do.call(rbind,complete.records.list)
    
    # Just to be safe, sort by month then year
    complete.records.df <- complete.records.list[order(complete.records.list$YEAR,complete.records.list$MONTH),]
    
    # Split by gaps in years
    # We want continuous sequences of years for linear interpolation
    splits <- cumsum(c(F,diff(complete.records.df$YEAR)>1))
    complete.records.df.split <- split(complete.records.df,splits)
    
    # (3) linearly interpolate across gaps <= the cutoff length
    interp.records.df.split <- lapply(complete.records.df.split,function(data.split){
      # Get matrix of only data
      data.matrix <- as.matrix(data.split[,4:34])
      
      # Get the number of days per month in the records
      n.days <- monthDays(as.Date(paste(data.split$YEAR,data.split$MONTH,'01',sep='-')))
      
      ## Unnwrap each row, accounting for number of days in the month
      data.matrix.unwrapped <- unwrapRows(data.matrix,n.days)
      
      # Interpolate
      data.matrix.unwrapped.interp <- round(na.approx(data.matrix.unwrapped, na.rm = F))
      
      # coerce back into a matrix
      data.matrix.final <- rewrapRows(data.matrix.unwrapped.interp,n.days)
      
      # join with year/month and return
      return(cbind(data.split[,c('YEAR','MONTH')],data.matrix.final))
    })
    
    # bind all chuncks
    interp.records.final <- do.call(rbind,interp.records.df.split)
    
    ## (4) fill by average those days at the beginning and end of years with missing data
    # split by year
    interp.records.final.years <- split(interp.records.final,interp.records.final$YEAR)
    # stack in 3D array
    interp.records.final.years.stack <- do.call(abind,c(interp.records.final.years,list(along=0))) # Gives 2 x 4 x 5
    class(interp.records.final.years.stack) <- "integer"
    # Calculate the mean of all years
    interp.records.final.years.stack.mean <- apply(interp.records.final.years.stack,c(2,3),function(X){mean(X,na.rm=T)})
    # Replace NA values with the mean of all years
    interp.records.final.years <- lapply(interp.records.final.years,function(year){
      year[is.na(year)] <- interp.records.final.years.stack.mean[is.na(year)]
      return(year)
    })
    # rbind
    interp.records.final <- do.call(rbind,interp.records.final.years)
    
    # And return
    return(interp.records.final)
  })
  
  if(any(sapply(data.list, is.null))){
    return(NULL)
  }
  
  ## Make sure you have the same years/months for each matrix
  yearMonth <- lapply(data.list,function(x){
    paste('Y',x$YEAR,'M',x$MONTH,sep='')
  })
  
  yearMonth.common <- Reduce(base::intersect,yearMonth)
  
  data.list <- lapply(data.list,function(x){
    this.yearMonth <- paste('Y',x$YEAR,'M',x$MONTH,sep='')
    return(x[match(yearMonth.common,this.yearMonth),])
  })
  
  ## Finally, (4) remove stations with less than 10 years of complete data (like WorldClim)
  if(is.null(data.list[[1]]) | nrow(data.list[[1]]) < (min.years*12)){
    return(NULL)
  }else{
    return(data.list)
  }  
}
