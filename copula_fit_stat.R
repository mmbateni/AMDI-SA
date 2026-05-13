# copula_fit_stat.R
library(copula)
library(ks)
source("empcop.R")
source("empkend.R")
source("copulafit2.R")

# ---- Settings ---------------------------------------------------------------
n_scale <- 12
overlap  <- n_scale - 1

mm <- dim(scor_spi_np)[3]
nn <- (nrow(scor_spi_np) * ncol(scor_spi_np)) - overlap

families <- c("Gaussian", "t", "Clayton", "Frank", "Gumbel", "Empirical")

# BUG FIX 3a: Classification breaks updated to 11-class mapping
y_edges_p <- c(0, 0.023, 0.055, 0.097, 0.212, 0.309, 0.691, 0.788, 0.903, 0.945, 0.977, 1)
y_edges   <- qnorm(y_edges_p)          
n_classes <- length(y_edges) - 1
mu_class  <- n_classes / 2              

# ---- Pre-allocate output arrays ---------------------------------------------
y_archem_K      <- array(NA, dim = c(nn, 3, mm))
K_c_theo        <- array(NA, dim = c(nn, 3, mm))
K_c_theo_KS     <- array(NA, dim = c(nn, 5, mm))
K_c_theo_MSE    <- array(NA, dim = c(nn, 3, mm))
best_cop_indx   <- matrix(NA,  nrow = nn, ncol = mm)
K_c_best        <- matrix(0,   nrow = nn, ncol = mm)  
index_cat       <- NULL
Rho_K_c_all     <- matrix(NA,  nrow = mm, ncol = 3)
cop_empri       <- matrix(NA,  nrow = nn, ncol = mm)
best_cop_all_KC <- vector("list", mm)
rmse_all        <- matrix(NA,  nrow = mm, ncol = 3)
Rho_copula_all  <- matrix(NA,  nrow = nn, ncol = 5)

# Helpers for theoretical Kendall functions
K_c_clayton_fn <- function(alpha, Ce) Ce * (1 + alpha - Ce^alpha) / alpha
K_c_frank_fn <- function(alpha, Ce) {
  if (abs(alpha) < 1e-8) return(Ce^2)
  Ce + ((1 - exp(-alpha * Ce)) / (alpha * exp(-alpha * Ce))) *
    log((1 - exp(-alpha)) / (1 - exp(-alpha * Ce)))
}
K_c_gumbel_fn <- function(alpha, Ce) Ce - (Ce * log(Ce)) / (alpha + 1)

# ---- Main loop --------------------------------------------------------------
for (ui in 1:mm) {
  
  percip_ind2 <- t(scor_spi_np[, , ui])
  u2 <- as.vector(percip_ind2)[n_scale:length(percip_ind2)]
  u1 <- u2[1:(length(u2) - lag)]
  
  ranking_SP1 <- rank(u1, ties.method = "average")
  AS_SP <- ranking_SP1 / (length(u1) + 1)
  
  smoist_ind2 <- t(score_smi_np[, , ui])
  v2 <- as.vector(smoist_ind2)[n_scale:length(smoist_ind2)]
  v1 <- v2[(lag + 1):length(v2)]
  
  ranking_SM1 <- rank(v1, ties.method = "average")
  AS_SM <- ranking_SM1 / (length(v1) + 1)
  
  w <- cbind(AS_SP, AS_SM)
  
  # ---- Fit copulas (MLE) ----------------------------------------------------
  fit_gauss   <- fitCopula(normalCopula(dim = 2),          data = w)
  fit_t       <- fitCopula(tCopula(dim = 2, df.fixed = FALSE), data = w)
  fit_clayton <- fitCopula(claytonCopula(dim = 2),         data = w)
  fit_frank   <- fitCopula(frankCopula(dim = 2),           data = w)
  fit_gumbel  <- fitCopula(gumbelCopula(dim = 2),          data = w)
  
  rho_gauss     <- coef(fit_gauss)
  params_t      <- coef(fit_t)
  rhohat        <- params_t["rho.1"]
  nuhat         <- params_t["df"]
  alpha_clayton <- coef(fit_clayton)
  alpha_frank   <- coef(fit_frank)
  alpha_gumbel  <- coef(fit_gumbel)
  
  # Theoretical copula CDFs at data points
  y_gauss   <- pCopula(w, fit_gauss@copula)
  y_t       <- pCopula(w, fit_t@copula)
  y_clayton <- pCopula(w, fit_clayton@copula)
  y_frank   <- pCopula(w, fit_frank@copula)
  y_gumbel  <- pCopula(w, fit_gumbel@copula)
  
  cop_empri[, ui] <- empcop(w)
  Ce <- cop_empri[, ui]
  
  # BUG FIX 3d: Call empkend on Gaussian and t CDFs instead of passthrough
  K_c_gauss_ks   <- empkend(w, y_gauss)
  K_c_t_ks       <- empkend(w, y_t)
  K_c_clay_ks    <- K_c_clayton_fn(alpha_clayton, Ce)
  K_c_frnk_ks    <- K_c_frank_fn(alpha_frank,   Ce)
  K_c_gumb_ks    <- K_c_gumbel_fn(alpha_gumbel, Ce)
  
  K_c_theo_KS[, , ui] <- cbind(K_c_gauss_ks, K_c_t_ks,
                               K_c_clay_ks,  K_c_frnk_ks, K_c_gumb_ks)

  # BUG FIX 3b & 3c: AIC Selection logic matched to MATLAB
  aic_ui <- c(AIC(fit_gauss), AIC(fit_t), AIC(fit_clayton), AIC(fit_frank), AIC(fit_gumbel))
  I_sorted <- order(aic_ui)
  I_cand_best <- I_sorted[1] 
  I_cand_2nd  <- I_sorted[2] 
  
  K_c_best[, ui]        <- K_c_theo_KS[, I_cand_best, ui]
  best_cop_all_KC[[ui]] <- families[I_cand_best]
  
  # ---- Transform to normal scores and discretise ----------------------------
  input_Cop          <- K_c_best[, ui]
  V                  <- qnorm(input_Cop)
  
  # BUG FIX 3e: Clip infinities to match MATLAB
  V[V == -Inf] <- -4
  V[V == Inf]  <- 4
  best_cop_indx[, ui] <- V
  
  index_cat_1 <- cut(best_cop_indx[, ui], breaks = y_edges, include.lowest = TRUE)
  index_cat   <- cbind(index_cat, index_cat_1)
}

index_cat <- sweep(index_cat, 2, mu_class, `-`)

save(y_archem_K, K_c_theo, K_c_theo_KS, K_c_theo_MSE,
     best_cop_indx, index_cat, Rho_K_c_all,
     cop_empri, best_cop_all_KC,
     rmse_all, Rho_copula_all,
     file = "coupla_derived.Rdata")
