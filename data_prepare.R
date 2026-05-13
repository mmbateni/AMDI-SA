# data_prepare.R
# Transforms 3-D SMI score array and computes autocorrelation for each basin.

# ---- Flatten and trim 3-D SMI scores ----------------------------------------
bs_score_smi_np_trans <- matrix(NA, nrow = 325, ncol = 79)

for (j in 1:79) {
  tmp <- as.vector(t(bs_score_smi_np[, , j]))
  
  # BUG FIX: MATLAB 12:end starts at element 12. 
  # Changed from 13:length(tmp) to 12:length(tmp)
  tmp <- tmp[12:length(tmp)]                      
  bs_score_smi_np_trans[, j] <- matrix(tmp, ncol = 1)
}

# ---- Autocorrelation for each basin -----------------------------------------
n_obs         <- nrow(bs_score_smi_np_trans)  
acf_matrix    <- matrix(NA, nrow = 12, ncol = 79)
lags_matrix   <- matrix(NA, nrow = 12, ncol = 79)
bounds_matrix <- matrix(NA, nrow = 12, ncol = 79)

for (i in 1:79) {
  result <- acf(bs_score_smi_np_trans[, i], lag.max = 12, plot = FALSE)
  
  # Drop lag-0 (index 1, always = 1.0) and keep lags 1:12.
  acf_matrix[, i]  <- result$acf[2:13, 1, 1]
  lags_matrix[, i] <- result$lag[2:13, 1, 1]
  
  # Approximate 95% confidence bounds
  bounds_matrix[, i] <- rep(qnorm(0.975) / sqrt(n_obs), 12)
}
