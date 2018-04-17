## A function that rewraps a vector into a matrix
## given a vector (n)  specifying the length of rows.
## n can be either a constant (in which case it will be repeated), or a vector
rewrapRows <- function(vec,n){
  list.out <- vector('list',length(n))
  n.sums <- c(0,cumsum(n))
  for(i in 2:length(n.sums)){
    temp.out <- vector('numeric',31)
    temp.out[] <- NA
    temp.out[1:length((n.sums[i-1]+1):n.sums[i])] <- vec[(n.sums[i-1]+1):n.sums[i]]
    list.out[[i]] <- temp.out
  }
  return(do.call(rbind,list.out))
}