copulafit2 <- function(u) {
  require(pracma)
  d <- ncol(u)
  
  if (ndim(u) != 2 || d < 2) {
    stop("InvalidDataDimensions")
  } else if (!all(u > 0 & u < 1)) {
    stop("DataOutOfRange")
  }
  
  families <- c('clayton', 'frank', 'gumbel')
  alphaHat <- numeric(length(families))
  rmse <- numeric(length(families))
  
  cop_empri <- empcop(u)
  K_c <- empkend(u, cop_empri)
  
  for (fa in seq_along(families)) {
    FAMILY <- families[fa]
    family <- match.arg(FAMILY, families)
    
    switch(family,
           clayton = {
             nloglf <- function(alpha) rmse_clayton(alpha, cop_empri, K_c)
             lowerBnd <- options$TolBnd
           },
           frank = {
             nloglf <- function(alpha) rmse_frank(alpha, cop_empri, K_c)
             lowerBnd <- bracket1D(nloglf, 5, -5)
             if (!is.finite(lowerBnd)) {
               stop("NoLowerBnd")
             }
           },
           gumbel = {
             nloglf <- function(alpha) rmse_gumbel(alpha, cop_empri, K_c)
             lowerBnd <- 1 + options$TolBnd
           })
    
    c(lowerBnd, upperBnd) <- bracket1D(nloglf, lowerBnd, 5)
    if (!is.finite(upperBnd)) {
      stop("NoUpperBnd")
    }
    
    alphaHat[fa] <- optimize(nloglf, c(lowerBnd, upperBnd), maximum = FALSE)$minimum
    rmse[fa] <- c(rmse_clayton(alphaHat[1], cop_empri, K_c),
                  rmse_frank(alphaHat[2], cop_empri, K_c),
                  rmse_gumbel(alphaHat[3], cop_empri, K_c))
  }
  
  return(list(alphaHat = alphaHat, rmse = rmse))
}
####################
