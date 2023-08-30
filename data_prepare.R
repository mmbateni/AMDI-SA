library(astsa) # Load the 'astsa' package for the 'autocorr' function

bs_score_smi_np_trans <- matrix(NA, nrow = 325, ncol = 79)
for (j in 1:79) {
  temp_bs_score_smi_np <- bs_score_smi_np[, , j]
  temp_bs_score_smi_np <- t(temp_bs_score_smi_np)
  temp_bs_score_smi_np <- as.vector(temp_bs_score_smi_np)
  temp_bs_score_smi_np <- temp_bs_score_smi_np[13:length(temp_bs_score_smi_np)]
  temp_bs_score_smi_np <- matrix(temp_bs_score_smi_np, ncol = 1)
  bs_score_smi_np_trans[, j] <- temp_bs_score_smi_np
}

acf_matrix <- matrix(NA, nrow = 12, ncol = 79)
lags_matrix <- matrix(NA, nrow = 12, ncol = 79)
bounds_matrix <- matrix(NA, nrow = 12, ncol = 79)

for (i in 1:79) {
  result <- autocorr(bs_score_smi_np_trans[, i], lag.max = 12)
  acf_matrix[, i] <- result$acf
  lags_matrix[, i] <- result$lag
  bounds_matrix[, i] <- result$ci
}
