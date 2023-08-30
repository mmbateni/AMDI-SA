y_edges2_p <- c(0, 0.05, 0.2, 0.3, 1)
# For 4 classes: c(0, 0.05, 0.2, 0.3, 1)
# Corresponding quantiles for 7 classes: c(0.000, 0.023, 0.067, 0.159, 0.841, 0.933, 0.977, 1.000)

# Convert quantiles to normal space
y_edges2 <- qnorm(y_edges2_p, mean = 0, sd = 1)

# Discretize using the quantiles
index_cat_mc <- cut(best_cop_indx, breaks = y_edges2, right = TRUE)

muu <- length(y_edges2) / 2
index_cat <- (-1) * (as.numeric(index_cat_mc) - muu) + muu
# Now, 1 means wettest and 4 (end) means driest

# Save index_cat to a file
saveRDS(index_cat, file = "coupla_derived.RDS", append = TRUE)
