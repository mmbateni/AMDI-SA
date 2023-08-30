empkend <- function(x, tq) {
  
  n <- nrow(x)
  n_tq <- length(tq)
  
  w <- numeric(n)
  wu <- matrix(NA, nrow=n, ncol=2)
  wu[,1] <- sort(x[,1])
  wu[,2] <- sort(x[,2])
  
  hi_ind <- matrix(NA, nrow=n, ncol=1) 
  hj_ind <- matrix(NA, nrow=n, ncol=1)
  
  Kn <- numeric(n_tq)
  
  for(k in 1:n_tq) {
    for(j in 1:n) {
      for(i in 1:n) {
        hi_ind[i] <- (x[i,1] < wu[j,1]) & (x[i,2] < wu[j,2]) 
      }
      w[j] <- sum(hi_ind)/(n+1)
      hj_ind[j] <- (w[j] < tq[k]) 
    }
    Kn[k] <- sum(hj_ind)/n
  }
  
  return(Kn)
  
}