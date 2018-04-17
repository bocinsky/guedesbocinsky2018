## A function that unwraps a matrix and only keeps the first n elements
## n can be either a constant (in which case it will be repeated), or a vector
unwrapRows <- function(mat,n){
  n <- rep_len(n,nrow(mat))
  i <- 0
  out <- lapply(1:nrow(mat),function(i){
    return(mat[i,1:n[i]])
  })
  return(as.numeric(do.call(c,out)))
}
