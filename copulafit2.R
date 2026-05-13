copulafit2 <- function(u) {
  require(pracma)
  
  TolBnd <- 1e-6   
  
  d <- ncol(u)
  if (!is.matrix(u) || d < 2) {
    stop("InvalidDataDimensions: u must be a matrix with at least 2 columns.")
  }
  if (!all(u > 0 & u < 1)) {
    stop("DataOutOfRange: all values of u must lie strictly in (0, 1).")
  }
  
  families  <- c("clayton", "frank", "gumbel")
  alphaHat  <- numeric(length(families))
  rmse_vals <- numeric(length(families))
  
  cop_empri <- empcop(u)
  K_c       <- empkend(u, cop_empri)
  
  # ---- RMSE helper functions ------------------------------------------------
  K_c_clayton <- function(alpha, Ce) {
    Ce * (1 + alpha - Ce^alpha) / alpha
  }
  K_c_frank <- function(alpha, Ce) {
    Ce + ((1 - exp(-alpha * Ce)) / (alpha * exp(-alpha * Ce))) *
      log((1 - exp(-alpha)) / (1 - exp(-alpha * Ce)))
  }
  K_c_gumbel <- function(alpha, Ce) {
    Ce - (Ce * log(Ce)) / (alpha + 1)
  }
  
  # BUG FIX: Objective function changed to match MATLAB's MAE-like approach
  # nll2 = sqrt((K_c(:) - K_c_pred(:)).^2); nll1 = sum(nll2) ./ sqrt(numel(K_c));
  rmse_fn <- function(Khat, Kemp) {
    sum(abs(Khat - Kemp)) / sqrt(length(Khat))
  }
  
  rmse_clayton <- function(alpha) rmse_fn(K_c_clayton(alpha, cop_empri), K_c)
  rmse_frank   <- function(alpha) rmse_fn(K_c_frank(alpha, cop_empri),   K_c)
  rmse_gumbel  <- function(alpha) rmse_fn(K_c_gumbel(alpha, cop_empri),  K_c)
  
  # ---- Optimise each family -------------------------------------------------
  for (fa in seq_along(families)) {
    family <- families[fa]
    
    if (family == "clayton") {
      nloglf   <- rmse_clayton
      lowerBnd <- TolBnd
      upperBnd <- 20          
      
    } else if (family == "frank") {
      nloglf   <- rmse_frank
      lowerBnd <- -20
      upperBnd <- 20
      
    } else if (family == "gumbel") {
      nloglf   <- rmse_gumbel
      lowerBnd <- 1 + TolBnd  
      upperBnd <- 20
    }
    
    res          <- optimize(nloglf, interval = c(lowerBnd, upperBnd))
    alphaHat[fa] <- res$minimum
    rmse_vals[fa] <- res$objective
  }
  
  return(list(
    families = families,
    alphaHat = alphaHat,
    rmse     = rmse_vals
  ))
}
