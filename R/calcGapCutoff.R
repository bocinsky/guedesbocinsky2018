## This function calculates the empirical cumulative distribution function,  
## then uses a user-defined probability level to calculate the "normal" gap length.
## That way, GHCN records (years) with unusually long gap lengths may be discarded.
calcGapCutoff <- function(rleVector, pLevel=0.95){
  f.all.rles <- ecdf(rleVector)
  x <- seq(0,31,length.out=100000)
  out <- mean(x[which(abs(f.all.rles(x)-pLevel)==min(abs(f.all.rles(x)-pLevel)))])
  cutoff <- ceil(out)
  return(cutoff)
}
