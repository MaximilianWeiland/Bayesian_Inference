# function for exporting a table as tex file
export_table <- function(df, destination) {
  posterior_table <- print(xtable(df),
                           type = "latex", include.rownames = FALSE,
                           print.results = FALSE)
  writeLines(posterior_table, destination)
}

# function for exporting a plot as pdf file
export_plot <- function(plot, destination) {
  pdf(destination, width = 8, height = 6)
  print(plot)
  dev.off()
}

# function to calculate waic based on posterior matrix
calculate_waic <- function(posterior) {
  
  # extract the posterior matrix for log likelihoods
  loglik_mat <- posterior$draws(format="matrix", variables="log_lik")
  
  # calculate lppd and p_waic
  lppd <- sum(log(apply(exp(loglik_mat),2,mean)))
  p_waic <- sum(apply(loglik_mat, 2, var))
  
  # based on these calculate the waic value and return
  WAIC <- -2 * (lppd - p_waic)
  WAIC
}

# function to check the performance of mh samples
check_performance_mh <- function(samples) {
  acceptance_rate <- mean(samples$accepted)
  cat(sprintf("Acceptance rate: %.2f%%\n", 100 * acceptance_rate))
}


# function to run and bind multiple mh chains
run_bind_mh_logreg <- function(n_chains = 3, X, y, n_samples = 1500,
                               warmup = 500, proposal_sd = 0.1) {
  
  # set up vector to store the chains in
  mh_chains <- vector("list", n_chains)
  
  # loop over chains
  for (i in 1:n_chains) {
    
    # set seed and produce different starting values
    set.seed(100 + i)
    start_beta <- rnorm(ncol(X), 0, 1)
    
    # run the chain and store result
    mh_chains[[i]] <- mh(
      X, y, n_samples = n_samples, warmup = warmup, proposal_sd = proposal_sd,
      start_beta = start_beta
    )
  }
  
  # stack all parameters
  beta_combined <- do.call(rbind, lapply(mh_chains, function(out) out$beta))
  ll_combined <- unlist(lapply(mh_chains, function(out) out$ll))
  yrep_combined <- do.call(rbind, lapply(mh_chains, function(out) out$y_rep))
  accepted_combined <- unlist(lapply(mh_chains, function(out) out$accepted))
  
  # combine to one compound list and return
  post_list <- list(
    beta = beta_combined,
    ll = ll_combined,
    y_rep = yrep_combined,
    accepted = accepted_combined
  )
  post_list
}

# extended version to also handle results of the multilevel model
bind_mh_multilevel <- function(n_chains = 3, X, U, y, country, n_samples = 1500,
                               warmup = 500, proposal_sd = 0.1) {
  
  # set up vector to store the chains in
  mh_chains <- vector("list", n_chains)
  
  # loop over chains
  for (i in 1:n_chains) {
    
    # run the chain and store result
    mh_chains[[i]] <- mh(X, U, y, country, n_samples,
                         warmup, proposal_sd)
  }
  
  # stack all parameters
  beta_combined <- abind::abind(lapply(mh_chains, function(out) out$beta), along = 1)
  gamma_combined <- abind::abind(lapply(mh_chains, function(out) out$gamma), along = 1)
  sigma_beta_combined <- do.call(rbind, lapply(mh_chains, function(out) out$sigma_beta))
  ll_combined <- unlist(lapply(mh_chains, function(out) out$ll))
  yrep_combined <- do.call(rbind, lapply(mh_chains, function(out) out$y_rep))
  accepted_combined <- unlist(lapply(mh_chains, function(out) out$accepted))
  
  # combine to one compound list
  post_list <- list(
    beta = beta_combined,
    gamma = gamma_combined,
    sigma_beta = sigma_beta_combined,
    ll = ll_combined,
    y_rep = yrep_combined,
    accepted = accepted_combined
  )
  post_list
}

# function to visualize traceplots
traceplots <- function(posterior_mat, n_chains = 3, iterations = 1000) {
  
  # define dimensionality variables
  chain_length <- iterations
  params <- colnames(posterior_mat)
  n_params <- length(params)
  
  # set the colors
  colors <- c("deepskyblue", "dodgerblue3", "steelblue")
  
  # set the layout
  par(mfrow = c(1, n_params), mar = c(4, 4, 2, 1))
  
  # split chains according to the chain length (3 x 1000)
  chain_list <- lapply(1:n_chains, function(i) {
    start <- (i - 1) * chain_length + 1
    end <- i * chain_length
    posterior_mat[start:end, , drop = FALSE]
  })
  
  # plot all parameters for all chains
  for (j in 1:n_params) {
    
    # plot first chain and define plot settings
    plot(chain_list[[1]][, j], type = "l", col = colors[1],
         ylab = "", xlab = "", main = params[j],
         ylim = range(posterior_mat[, j]))
    
    # loop over remaining chains and plot them with respective color
    for (i in 2:n_chains) {
      lines(chain_list[[i]][, j], col = colors[i])
    }
  }
}

# function for extracting all parameter draws from hmc posterior matrix
extract_matrix_hmc <- function(posterior, which_columns, column_names) {
  
  # get the posterior matrix and set new column names
  post_mat <- posterior$draws(format="matrix")[, which_columns]
  colnames(post_mat) <- column_names
  post_mat
}

# function for extracting all parameter draws from mh posterior list
extract_matrix_mh <- function(posterior_list, column_names) {
  
  # get elements from the list
  ll <- posterior_list$ll
  beta <- posterior_list$beta
  post_mat <- cbind(ll, beta)
  
  # set new column names to posterior matrix
  colnames(post_mat) <- column_names
  post_mat
}

# function to calculate population average effects
pop_avg_effects <- function(X, posterior, var) {
  
  # get the index of the variable for which marginal effects should be plotted
  var_index <- which(colnames(X) == var)
  
  # define which values to consider with special handling of the dummy variable
  if (length(unique(X[, var_index])) == 2 && all(unique(X[, var_index]) %in% c(0, 1))) {
    # dummy variable only gets 0 and 1
    grid_vals <- c(0, 1)
  } else {
    grid_vals <- seq(-2, 2, length.out = 50)
  }
  
  # set up a df in which to store the results
  n_iter <- nrow(posterior)
  n_obs <- nrow(X)
  pred_summary <- data.frame(value = grid_vals, mean = NA, lower = NA, upper = NA)
  
  # loop over all x-values of the variable of interest
  for (i in seq_along(grid_vals)) {
    
    # create a copy of the data and manipulate the variable of interest
    X_tmp <- X
    X_tmp[, var_index] <- grid_vals[i]
    
    # predict values for each mcmc sample and store p-value in matrix
    pred_probs <- matrix(NA, nrow = n_iter, ncol = n_obs)
    for (j in 1:n_iter) {
      pred_probs[j, ] <- 1 / (1 + exp(-X_tmp %*% matrix(posterior[j, ], ncol = 1)))
    }
    
    # take the average prediction for each sample iteration
    avg_probs <- rowMeans(pred_probs)
    
    # calculate summary statistics across sample iterations
    pred_summary$mean[i] <- mean(avg_probs)
    pred_summary$lower[i] <- quantile(avg_probs, 0.025)
    pred_summary$upper[i] <- quantile(avg_probs, 0.975)
  }
  pred_summary
}

# function to create a long summary of parameter values
create_long_summary <- function(post_matrix, parameter, country_names, var_names) {
  
  # get the beta draws
  par_draws <- post_matrix[,grep(parameter, colnames(post_matrix)), drop = FALSE]
  
  # set up a summary df
  summary_df <- data.frame(
    country = integer(length = ncol(post_matrix)),
    variable = integer(length = ncol(post_matrix)),
    mean = numeric(ncol(post_matrix)),
    lower = numeric(ncol(post_matrix)),
    upper = numeric(ncol(post_matrix)),
    stringsAsFactors = FALSE
  )
  
  # get the indices stored in the column names of the matrix
  colnames_split <- regmatches(colnames(post_matrix),
                               gregexpr("[0-9]+", colnames(post_matrix)))
  
  # explicitly store row and column indices
  row_indices <- as.integer(sapply(colnames_split, function(x) x[1]))
  col_indices <- as.integer(sapply(colnames_split, function(x) x[2]))
  
  # get the correct country and variable names based on indices and store in df
  summary_df$country <- country_names[row_indices]
  summary_df$variable <- var_names[col_indices]
  
  # move over posterior matrix columns and calculate summary statistics
  for (i in 1:ncol(post_matrix)) {
    samples <- post_matrix[, i]
    summary_df$mean[i]  <- mean(samples)
    summary_df$lower[i] <- quantile(samples, 0.025)
    summary_df$upper[i] <- quantile(samples, 0.975)
  }
  summary_df
}

# the same function, but for mh-posterior
create_long_summary_mh <- function(post_array, country_names, var_names) {
  
  # get dimensions of the array
  n_samples <- dim(post_array)[1]
  n_countries <- dim(post_array)[2]
  n_vars <- dim(post_array)[3]
  
  # create the summary df
  summary_df <- data.frame(
    country = rep(country_names, each = n_vars),
    variable = rep(var_names, times = n_countries),
    mean = NA,
    lower = NA,
    upper = NA,
    stringsAsFactors = FALSE
  )
  
  # loop over all countries and all variables and calculate summary
  idx <- 1
  for (j in 1:n_countries) {
    for (k in 1:n_vars) {
      draws <- post_array[, j, k]
      summary_df$mean[idx]  <- mean(draws)
      summary_df$lower[idx] <- quantile(draws, 0.025)
      summary_df$upper[idx] <- quantile(draws, 0.975)
      idx <- idx + 1
    }
  }
  summary_df
}

# function to plot the point estimates with confidence intervals
plot_estimates_ci <- function(posterior_summary, country_names, individual_predictors) {
  
  # get dimensionality of country names, number of models and the type of models
  n_countries <- length(country_names)
  model_levels <- unique(posterior_summary$model_type)
  n_models <- length(model_levels)
  
  # define colors for plotting and model names
  model_colors <- c("black", "red")
  names(model_colors) <- model_levels
  model_pchs <- c(19, 17)
  names(model_pchs) <- model_levels
  
  # set the plot layout
  par(mfrow = c(1, length(individual_predictors)), mar = c(4, 6, 3, 1), oma = c(4, 0, 2, 0))
  
  # loop over all predictors and create a temporary subset of the df
  for (pred in individual_predictors) {
    sub <- posterior_summary[posterior_summary$variable == pred, ]
    
    # set country to a factor variable and order according to country
    sub$country <- factor(sub$country, levels = rev(country_names))
    sub <- sub[order(sub$country, decreasing = TRUE), ]
    
    # get country and model index
    country_idx <- as.integer(sub$country)
    model_idx <- match(sub$model_type, model_levels)
    
    # set the y-position for each country for each model type
    y_pos <- country_idx * n_models - (model_idx - 1)
    
    # define x and y limit ranges
    ylim_range <- c(0.5, n_countries * n_models + 0.5)
    xlim_range <- range(sub$lower, sub$upper)
    
    # plot the point estimate of each country's slope estimate
    plot(x = sub$mean, y = y_pos,
         xlim = xlim_range,
         ylim = ylim_range,
         xlab = "Estimate",
         ylab = "",
         yaxt = "n",
         main = pred,
         pch = model_pchs[sub$model_type],
         col = model_colors[sub$model_type])
    
    # add the confidence intervals
    segments(x0 = sub$lower, x1 = sub$upper, y0 = y_pos, y1 = y_pos,
             col = model_colors[sub$model_type], lwd = 2)
    
    # add axis ticks and labels (the country names)
    axis_ticks <- seq(n_models / 2 + 0.5, n_countries * n_models, by = n_models)
    axis_labels <- rev(country_names)
    axis(2, at = axis_ticks, labels = axis_labels, las = 2)
    
    # add vertical dashed line at zero only if 0 in x-axis range
    if (xlim_range[1] <= 0 && xlim_range[2] >= 0) {
      segments(x0 = 0, y0 = ylim_range[1], x1 = 0, y1 = ylim_range[2], 
               lty = 2, col = "gray")
    }
  }
}

# function to create a dense summary for the gamma values
create_gamma_summary <- function(gamma_mat, group_var_names, indiv_var_names) {
  
  # create an empty matrix in which gamma summaries will be stored
  result_mat <- matrix(NA, nrow = length(group_var_names), ncol = length(indiv_var_names),
                       dimnames = list(group_var_names, indiv_var_names))

  # extract column names and save the indices defined in them
  col_names <- colnames(gamma_mat)
  indices <- regmatches(col_names, gregexpr("[0-9]+", col_names))
  
  # loop over all column names
  for (idx in seq_along(col_names)) {
    
    # get all draws from the current column
    draws <- gamma_mat[, idx]
    
    # get group and individual level indices from the column names
    group_idx <- as.integer(indices[[idx]][1])
    individual_idx <- as.integer(indices[[idx]][2])
    
    # calculate summary stats
    mean <- mean(draws)
    sd <- sd(draws)
    
    # store summary in the correct cell
    result_mat[group_idx, individual_idx] <- sprintf("%.2f (%.2f)", mean, sd)
  }
  
  # convert the result to a dataframe
  result_df <- as.data.frame(result_mat, stringsAsFactors = FALSE)
  
  # get the group variable as separate variable
  result_df <- cbind(`Group Variable` = rownames(result_df), result_df)
  result_df
}
