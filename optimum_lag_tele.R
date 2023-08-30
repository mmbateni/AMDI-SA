ind_xls <- read.xlsx("index_basin.xlsx", sheet = "ind", startRow = 2, endCol = 3)
ind_xls_t_bsn <- ind_xls[, 1]
n_bs <- length(ind_xls_t_bsn)
scale <- 1:12
n_scale <- length(scale)
overlap <- n_scale - 1
bc <- ncol(best_cop_indx)
bd <- nrow(best_cop_indx)

## Rest of the code involves mostly data manipulation and calculations.
## Please make sure to use appropriate functions in R for file reading and saving.
