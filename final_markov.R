library(openxlsx)

ind_xls <- read.xlsx("index_basin.xlsx", sheet = "ind", startRow = 2, endCol = 3)
ind_basin <- ind_xls[, 1]
n_c <- 4
n_ci <- 3
nm <- n_ci * n_c
mm <- nrow(index_cat)
nn <- ncol(index_cat)
m <- ceiling(mm * (2 / 3))
l_q <- (mm - m) + 1

## Rest of the code involves mostly data manipulation and calculations.
## Please make sure to use appropriate functions in R for file reading and saving.
