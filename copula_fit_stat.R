# Load required packages
library(copula)
library(ks)

# Initialize variables
n_scale <- 12
overlap <- n_scale - 1
mm <- size(scor_spi_np, 3)
nn <- (nrow(scor_spi_np) * ncol(scor_spi_np)) - overlap  

families <- c("Gaussian", "t", "Clayton", "Frank", "Gumbel", "Empirical")

# Preallocate output matrices
y_archem_K <- array(NA, dim=c(nn, 3, mm)) 
K_c_theo <- array(NA, dim=c(nn, 3, mm))
K_c_theo_KS <- array(NA, dim=c(nn, 5, mm)) 
K_c_theo_MSE <- array(NA, dim=c(nn, 3, mm))
best_cop_indx <- matrix(NA, nrow=nn, ncol=mm)
index_cat <- NULL 
Rho_K_c_all <- matrix(NA, nrow=mm, ncol=3)
cop_empri <- matrix(NA, nrow=nn, ncol=mm)
H_KC_ks <- matrix(NA, nrow=mm, ncol=5)
P_KC_ks <- matrix(NA, nrow=mm, ncol=5) 
best_cop_all_KC <- vector("list", mm)
rmse_all <- matrix(NA, nrow=mm, ncol=3)
H_KC_ks_MSE <- matrix(NA, nrow=mm, ncol=3) 
P_KC_ks_MSE <- matrix(NA, nrow=mm, ncol=3)
Rho_copula_all <- matrix(NA, nrow=nn, ncol=5)
best_cop_all <- vector("list", mm)
ind_H_KC <- matrix(NA, nrow=mm, ncol=1)

for(ui in 1:mm){
  
  # Extract data
  percip_ind1 <- scor_spi_np[, , ui]
  percip_ind2 <- t(percip_ind1)
  u2 <- percip_ind2[n_scale:length(percip_ind2)]
  u1 <- u2[1:(length(u2)-lag)]
  
  # Get ranks
  ranking_SP1 <- rank(u1, ties.method="average")
  AS_SP <- ranking_SP1 / (length(u1) + 1)
  
  smoist_ind1 <- score_smi_np[, , ui] 
  smoist_ind2 <- t(smoist_ind1)
  v2 <- smoist_ind2[n_scale:length(smoist_ind2)]
  v1 <- v2[lag+1:length(v2)]
  ranking_SM1 <- rank(v1, ties.method="average")
  AS_SM <- ranking_SM1 / (length(v1) + 1)
  
  w <- cbind(AS_SP, AS_SM)
  
  # Fit copulas
  fit_gauss <- fitCopula(normalCopula(dim = 2), data = w)
  Rho1 <- fit_gauss@rho
  Rho <- mean(diag(Rho1))
  fit_t <- fitCopula(tCopula(dim = 2, df = 4), data = w)
  rhohat1 <- fit_t@rho
  rhohat <- mean(diag(rhohat1))
  nuhat <- fit_t@df
  
  fit_clayton <- fitCopula(claytonCopula(dim = 2), data = w)
  alpha_clayton <- fit_clayton@theta
  
  fit_frank <- fitCopula(frankCopula(dim = 2), data = w)
  alpha_frank <- fit_frank@theta
  
  fit_gumbel <- fitCopula(gumbelCopula(dim = 2), data = w)
  alpha_gumbel <- fit_gumbel@theta
  
  y_gauss <- pcopula(fit_gauss, w)
  y_t <- pcopula(fit_t, w)
  y_clayton <- pcopula(fit_clayton, w) 
  y_frank <- pcopula(fit_frank, w)
  y_gumbel <- pcopula(fit_gumbel, w)
  
  # Fit by Kendall's tau
  fit_clayton_ks <- fitCopula(claytonCopula(dim = 2), data = w, method = "itau")
  fit_frank_ks <- fitCopula(frankCopula(dim = 2), data = w, method = "itau")
  fit_gumbel_ks <- fitCopula(gumbelCopula(dim = 2), data = w, method = "itau")
  
  Rho_K_c <- c(fit_clayton_ks@theta, fit_frank_ks@theta, fit_gumbel_ks@theta)
  
  # Compute empirical copula 
  u <- ecdf(w)(w)
  cop_empri[,ui] <- u[,1] * u[,2]
  
  # Compute Kendall's tau
  C_tau1 <- cor(w[,1], w[,2], method="kendall")
  
  # Compute theoretical Kendall's functions
  K_c_clayton <- y_clayton * ((1 + alpha_clayton - (y_clayton^alpha_clayton))/alpha_clayton)
  ex_K_frank <- exp(-alpha_frank*y_frank)
  K_c_frank <- y_frank + (((1 - ex_K_frank)/(alpha_frank*ex_K_frank)) * 
                            (log((1 - exp(-alpha_frank))/(1 - ex_K_frank))))
  K_c_gumbel <- y_gumbel - ((y_gumbel*log(y_gumbel))/(alpha_gumbel + 1))
  K_c_theo[, , ui] <- rbind(K_c_clayton, K_c_frank, K_c_gumbel)
  
  K_c_clayton_ks <- cop_empri[, ui] * ((1 + alpha_clayton - (cop_empri[, ui]^alpha_clayton))/alpha_clayton)
  K_c_frank_ks <- cop_empri[, ui] + (((1 - exp(-alpha_frank*cop_empri[, ui]))/(alpha_frank*exp(-alpha_frank*cop_empri[, ui]))) *
                                       (log((1 - exp(-alpha_frank))/(1 - exp(-alpha_frank*cop_empri[, ui]))))) 
  K_c_gumbel_ks <- cop_empri[, ui] - ((cop_empri[, ui]*log(cop_empri[, ui]))/(alpha_gumbel + 1))
  K_c_gauss_ks <- pnorm(qnorm(y_gauss)) # Gaussian is just normal CDF
  K_c_t_ks <- pt(qt(y_t, df=nuhat), df=nuhat) # t is just t CDF
  K_c_theo_KS[, , ui] <- rbind(K_c_gauss_ks, K_c_t_ks, K_c_clayton_ks, K_c_frank_ks, K_c_gumbel_ks)
  
  # Compute empirical Kendall's function
  K_c <- empkend(w) 
  
  # Compute theoretical KC (best fit)
  K_c_clayton_ks <- cop_empri[, ui] * ((1 + alpha_clayton - (cop_empri[, ui]^alpha_clayton))/alpha_clayton)
  K_c_frank_ks <- cop_empri[, ui] + (((1 - exp(-alpha_frank*cop_empri[, ui]))/(alpha_frank*exp(-alpha_frank*cop_empri[, ui]))) *
                                       (log((1 - exp(-alpha_frank))/(1 - exp(-alpha_frank*cop_empri[, ui])))))
  K_c_gumbel_ks <- cop_empri[, ui] - ((cop_empri[, ui]*log(cop_empri[, ui]))/(alpha_gumbel + 1)) 
  K_c_gauss_ks <- pnorm(qnorm(y_gauss))
  K_c_t_ks <- pt(qt(y_t, df=nuhat), df=nuhat)
  K_c_theo_KS[, , ui] <- rbind(K_c_gauss_ks, K_c_t_ks, K_c_clayton_ks, K_c_frank_ks, K_c_gumbel_ks)
  
  # Compute theoretical KC (MSE fit)
  K_c_clayton_mse <- cop_empri[, ui] * ((1 + Rho_K_c[1] - (cop_empri[, ui]^Rho_K_c[1]))/Rho_K_c[1])
  K_c_frank_mse <- cop_empri[, ui] + (((1 - exp(-Rho_K_c[2]*cop_empri[, ui]))/(Rho_K_c[2]*exp(-Rho_K_c[2]*cop_empri[, ui]))) * 
                                        (log((1 - exp(-Rho_K_c[2]))/(1 - exp(-Rho_K_c[2]*cop_empri[, ui])))))
  K_c_gumbel_mse <- cop_empri[, ui] - ((cop_empri[, ui]*log(cop_empri[, ui]))/(Rho_K_c[3] + 1))
  K_c_theo_MSE[, , ui] <- rbind(K_c_clayton_mse, K_c_frank_mse, K_c_gumbel_mse)
  
  # Compare KC functions
  for(ks in 1:5){
    test <- ks.test(K_c_theo_KS[,ks,ui], K_c[,ui])
    H_KC_ks[ui, ks] <- as.numeric(test$p.value < 0.05)
    P_KC_ks[ui, ks] <- test$p.value
  }
  
  for(kse in 1:3){
    test <- ks.test(K_c_theo_MSE[,kse,ui], K_c[,ui])  
    H_KC_ks_MSE[ui, kse] <- as.numeric(test$p.value < 0.05)
    P_KC_ks_MSE[ui, kse] <- test$p.value
  }
  H_KC <- !(H_KC_ks[ui,])
  
  I_min <- which.min(rmse_all[ui,])
  I_min_2 <- which(sort(rmse_all[ui,]) == 2)
  I_min_3 <- which(sort(rmse_all[ui,]) == 3)
  
  H_KC_MSE <- !(H_KC_ks_MSE[ui,])
  
  K_c_best[,ui] <- K_c[,ui]
  best_cop_all_KC[[ui]] <- families[[6]]
  
  if(H_KC[I_best]){
    if(I_best > 2){
      K_c_best[,ui] <- K_c_theo[, I_best-2, ui] 
      best_cop_all_KC[[ui]] <- families[[I_best]]
    } else {
      K_c_best[,ui] <- K_c_theo[, I_best, ui]
      best_cop_all_KC[[ui]] <- families[[I_best]] 
    }
  }
  
  if(all(K_c_best[,ui] == 0)){
    if(H_KC_MSE[I_min_3]){
      K_c_best[,ui] <- K_c_theo_MSE[, I_min_3, ui]
      best_cop_all_KC[[ui]] <- paste0(families[[I_min_3]], " LD")
    } else if(H_KC_MSE[I_min_2]){
      K_c_best[,ui] <- K_c_theo_MSE[, I_min_2, ui]
      best_cop_all_KC[[ui]] <- paste0(families[[I_min_2]], " LD")
    } else if(H_KC_MSE[I_min]){
      K_c_best[,ui] <- K_c_theo_MSE[, I_min, ui]
      best_cop_all_KC[[ui]] <- paste0(families[[I_min]], " LD")
    }
  }
  
  if(all(K_c_best[,ui] == 0)){
    K_c_best[,ui] <- K_c[,ui]
    best_cop_all_KC[[ui]] <- families[[6]]
  }
  # Transform to normal scores
  input_Cop <- K_c_best[,ui]
  best_cop_indx[,ui] <- qnorm(input_Cop)
  
  # Discretize
  index_cat_1 <- cut(best_cop_indx[,ui], breaks=y_edges, include.lowest=TRUE)
  index_cat <- cbind(index_cat, index_cat_1)
  
}

# Center discretized values
index_cat <- sweep(index_cat, 2, mu, `-`)

# Save outputs
save(y_archem_K, K_c_theo, K_c_theo_KS, K_c_theo_MSE, best_cop_indx, index_cat, Rho_K_c_all,  
     cop_empri, H_KC_ks, P_KC_ks, best_cop_all_KC, rmse_all, H_KC_ks_MSE, P_KC_ks_MSE,  
     Rho_copula_all, best_cop_all, ind_H_KC, file="coupla_derived.Rdata")