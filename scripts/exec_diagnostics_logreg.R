############### Run Non-Regularized Logistic Regression Models #################

# compile stan model
mod_non_reg <- cmdstan_model("model_files/logistic_regression.stan")

# setup the data list (valid for all models)
data_list <- list(N = length(y_germany),
                  K = ncol(X_ger_scaled),
                  X = X_ger_scaled,
                  y = y_germany
)

# provide a function that initializes random starting values for each chain
init_fun <- function() {
  list(
    beta = rnorm(ncol(X_ger_scaled), 0, 0.1)
  )
}

# run the sampler
post_logreg_nonreg_hmc <- mod_non_reg$sample(
  data = data_list,
  seed = 123,
  chains = 3,
  parallel_chains = 3,
  thin = 1,
  iter_warmup = 500,
  iter_sampling = 1000,
  init = init_fun,
  save_warmup = TRUE
)

# save the posterior matrix
post_matrix_logreg_nonreg_hmc <- post_logreg_nonreg_hmc$draws(format="matrix")
save(post_matrix_logreg_nonreg_hmc, file = "posteriors/hmc_logreg_noreg.RData")

# run the metropolis hastings algorithm for the non-regularized model

# load the model function
source("model_files/logistic_regression_mh.R")

# tune the standard deviation to get good results
proposal_sds_dummy <- c(0.2, 0.4, 0.6, 0.8, 1)

for (sd_dummy in proposal_sds_dummy) {
  # set up a vector of proposal sd = 0.1
  proposal_sd_vec <- rep(0.1, ncol(X_ger_scaled))
  proposal_sd_vec[4] <- sd_dummy
  
  result <- mh(X_ger_scaled, y_germany, n_samples = 1500,
               warmup = 500, proposal_sd_vec = proposal_sd_vec,
               start_beta = rep(1, ncol(X_ger_scaled)))
  
  cat("Chain with proposal sd for dummy variable:", sd_dummy, "\n")
  check_performance_mh(result)
  cat("----------------------------------------------\n")
}

# run three chains of the model with the tuned sd
optimal_sd_vec <- rep(0.1, ncol(X_ger_scaled))
optimal_sd_vec[4] <- 0.6
post_logreg_nonreg_mh <- run_bind_mh_logreg(n_chains = 3, X = X_ger_scaled,
                                         y = y_germany, n_samples = 1500,
                                         warmup = 500, proposal_sd = optimal_sd_vec)

# save the model
save(post_logreg_nonreg_mh, file = "posteriors/mh_logreg_noreg.RData")


################### Run Regularized Logistic Regression Models #################

# compile stan model
mod_reg <- cmdstan_model("model_files/logistic_regression_regularized.stan")

# run the sampler
post_logreg_reg_hmc <- mod_reg$sample(
  data = data_list,
  seed = 123,
  chains = 3,
  parallel_chains = 3,
  thin = 1,
  iter_warmup = 500,
  iter_sampling = 1000,
  init = init_fun,
  save_warmup = TRUE
)

# save the posterior matrix
post_matrix_logreg_reg_hmc <- post_logreg_reg_hmc$draws(format="matrix")
save(post_matrix_logreg_reg_hmc, file = "posteriors/hmc_logreg_reg.RData")

# run the metropolis hastings algorithm for the regularized model

# load the model function
source("model_files/logistic_regression_regularized_mh.R")

# tune the standard deviation to get good results
proposal_sds_dummy <- c(0.1, 0.2, 0.3, 0.4, 0.5)

for (sd_dummy in proposal_sds_dummy) {
  # set up a vector of proposal sd = 0.1
  proposal_sd_vec <- rep(0.1, ncol(X_ger_scaled))
  proposal_sd_vec[4] <- sd_dummy
  
  result <- mh(X_ger_scaled, y_germany, n_samples = 1500,
               warmup = 500, proposal_sd_vec = proposal_sd_vec,
               start_beta = rep(1, ncol(X_ger_scaled)))
  
  cat("Chain with proposal sd for dummy variable:", sd_dummy, "\n")
  check_performance_mh(result)
  cat("----------------------------------------------\n")
}

# run three chains of the model with the tuned sd
optimal_sd_vec <- rep(0.1, ncol(X_ger_scaled))
optimal_sd_vec[4] <- 0.4
post_logreg_reg_mh <- run_bind_mh_logreg(n_chains = 3, X = X_ger_scaled,
                                            y = y_germany, n_samples = 1500,
                                            warmup = 500, proposal_sd = optimal_sd_vec)
# save the model
save(post_logreg_reg_mh, file = "posteriors/mh_logreg_reg.RData")


######################## Check Convergence #####################################

# extract posterior matrices for all models and rename columns

# for the hmc models
column_names <- c("lp", "intercept", "satisfaction_gov",
                  "age", "unemployment")
post_mat_nonreg_hmc <- extract_matrix_hmc(post_logreg_nonreg_hmc, 1:5, column_names)
post_mat_reg_hmc <- extract_matrix_hmc(post_logreg_reg_hmc, 1:5, column_names)

# for the mh models
post_mat_nonreg_mh <- extract_matrix_mh(post_logreg_nonreg_mh, column_names)
post_mat_reg_mh <- extract_matrix_mh(post_logreg_reg_mh, column_names)

# arrange all models in a list and create vector for model names
model_list <- list(post_mat_nonreg_hmc, post_mat_nonreg_mh,
                   post_mat_reg_hmc, post_mat_reg_mh)
model_names <- c("HMC - Weakly Informative", "MH - Weakly Informative",
                 "HMC - Regularized", "MH - Regularized")

# loop over models, create traceplot and export
for (i in seq_along(model_list)) {
  pdf(paste0("figures/logreg_traceplot_", model_names[i], ".pdf"), width = 12, height = 3)
  traceplots(model_list[[i]], n_chains = 3, iterations = 1000)
  dev.off()
}


######################### Posterior Summaries ##################################

summary_list <- list()

for (i in seq_along(model_list)) {
  
  # extract the current posterior and remove ll
  posterior <- model_list[[i]]
  posterior <- posterior[, 2:ncol(posterior)]
  
  # calculate mean and sd
  means <- colMeans(posterior)
  sds <- apply(posterior, 2, sd)
  
  # format as mean (sd)
  mean_sd <- sprintf("%.2f (%.2f)", means, sds)
  
  # append to the list
  summary_list[[i]] <- c(Model = model_names[i], mean_sd)
}

# combine all entries into a df
summary_df <- do.call(rbind, summary_list)
summary_df <- as.data.frame(summary_df)

# set the column names
colnames(summary_df) <- c("Model", "intercept", "satisfaction_gov",
                          "age", "unemployment")

# export the summary table as a tex file
export_table(summary_df, "report/posterior_summary_logreg.tex")
