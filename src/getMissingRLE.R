## This function gets the run length encoding of joint NAs in the climate records.
## Useful for deciding the cutoff value for the number of contiguous days.
getMissingRLE <- function(vect){
  ## Read in both records
  vect.tmin <- vect[['TMIN']]
  vect.tmax <- vect[['TMAX']]
  
  ## Make sure they are both the same dimensions
  if(!equals(dim(vect.tmin),dim(vect.tmax))){
    return(NULL)
  }
  
  ## Get just the climate info
  vect.tmin.records <- as.matrix(vect.tmin[,3:33])
  vect.tmax.records <- as.matrix(vect.tmax[,3:33])
  
  ## Remove all records where we don't have both data
  vect.tmin.records[is.na(vect.tmax.records)] <- NA
  vect.tmax.records[is.na(vect.tmin.records)] <- NA
  
  ## Get the number of days per month in the records
  n.days <- monthDays(as.Date(paste(vect.tmin$YEAR,vect.tmin$MONTH,'01',sep='-')))
  
  ## Unnwrap each row, accounting for number of days in the month, and count
  ## consecutive missing values (these are the same)
  vect.tmin.records.rle <- rle(is.na(unwrapRows(vect.tmin.records,n.days)))
  #   vect.tmax.records.rle <- rle(is.na(unwrap.rows(vect.tmax.records,n.days)))
  
  out <- vect.tmin.records.rle$lengths[vect.tmin.records.rle$values]
  return(out)
}
