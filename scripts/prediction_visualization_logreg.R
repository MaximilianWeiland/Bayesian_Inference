# initialize list to store the predicted probabilities
pred_list <- list()

# list of model names for storing summaries under correct name
model_names <- c("nonreg", "reg")

# vector of variable names (names match the X-matrix)
variables <- c("satisfaction_gov_germany", "age_germany", "unemployment_germany")

# set up list with models for which predictions should be computed
model_list <- list(post_mat_nonreg_hmc, post_mat_reg_hmc)

# loop over models
for (m in seq_along(model_list)) {
  
  # extract the posterior and get only relevant variables
  posterior <- model_list[[m]][, 2:5]
  
  # loop over all variables in the posterior and calculate marginal effects
  for (var in variables) {
    pred_summary <- pop_avg_effects(X_ger_scaled, posterior, var)
    
    # store with descriptive name
    pred_list[[paste(model_names[m], var, sep = "_")]] <- pred_summary
  }
}

# specify variable names for plotting
var_names_plot <- c("satisfaction_gov", "age", "unemployment")

# logical vector in which is stored which variable is the dummy variable
dummy_bool <- c(FALSE, FALSE, TRUE)

# export the plot
pdf("figures/predictions_pop_average_logreg.pdf", width = 10, height = 5)

# set up grid for the plot
par(mfrow = c(1, 3))

# loop over all variables and build plot
for (i in seq_along(variables)) {
  
  # get the variable name and extract summary dfs
  variable <- variables[i]
  nonreg <- pred_list[[paste0("nonreg_", variable)]]
  reg <- pred_list[[paste0("reg_", variable)]]
  
  # get the variable name for plotting
  var_name <- var_names_plot[i]
  
  # check if it is the dummy variable
  dummy <- dummy_bool[i]
  
  # if it is not the dummy, plot the means as lineplot (first for non-regularized)
  if (dummy == F) {
    plot(
      nonreg$value, nonreg$mean, type = "l", col = "black", lwd = 2,
      ylim = range(c(nonreg$lower, nonreg$upper, reg$lower, reg$upper)),
      xlab = paste(var_name, "(z-standardized)"), ylab = "Predicted Probability",
      main = paste("Effect of", var_name)
    )
    
    # add confidence intervals
    polygon(
      c(nonreg$value, rev(nonreg$value)),
      c(nonreg$lower, rev(nonreg$upper)),
      col = adjustcolor("black", alpha.f = 0.2),
      border = NA
    )
    
    # add lines for the regularized model
    lines(reg$value, reg$mean, col = "red", lwd = 2)
    polygon(
      c(reg$value, rev(reg$value)),
      c(reg$lower, rev(reg$upper)),
      col = adjustcolor("red", alpha.f = 0.2),
      border = NA
    )
    
    # add a legend
    legend("topright", legend = c("Non-Regularized", "Regularized"),
           col = c("black", "red"), lwd = 2, bty = "n")
  }
  
  # plot for unemployment (which is the dummy variable)
  if (dummy == T) {
    
    # plot for the non-regularized model
    plot(
      x = nonreg$value, y = nonreg$mean, type = "p", pch = 16,
      ylim = range(c(nonreg$lower, nonreg$upper, reg$lower, reg$upper)),
      xlab = var_name, ylab = "Predicted Probability",
      xaxt = "n", main = paste("Effect of", var_name)
    )
    
    # add x-axis labels for dummy values
    axis(1, at = nonreg$value, labels = c("No", "Yes"))
    
    # add confidence intervals for non-regularized model
    arrows(nonreg$value, nonreg$lower, nonreg$value, nonreg$upper,
           angle = 90, code = 3, length = 0.05, col = "black")
    
    # add points for regularized model
    points(reg$value, reg$mean, col = "red", pch = 16)
    arrows(reg$value, reg$lower, reg$value, reg$upper,
           angle = 90, code = 3, length = 0.05, col = "red")
    
    # add a legend
    legend("topright", legend = c("Non-Regularized", "Regularized"),
           col = c("black", "red"), pch = 16, bty = "n")
  }
}
dev.off()
