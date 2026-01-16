############### Run Non-Regularized Logistic Regression Models #################

# compile stan model
mod_non_reg <- cmdstan_model("Task3/Stan_MH_Files/multilevel_logistic.stan")

# setup the data list (valid for all models)
data_list <- list(y = y_multilevel, 
                  J=J, 
                  N=N,
                  country=country,
                  X=X_multilevel_scaled,
                  U=U_scaled,
                  K=K,
                  L=L
)

# provide a function that initializes random starting values for each chain
init_fun <- function() {
  list(
    beta = matrix(rnorm(J * K, 0, 0.1), nrow = J),
    gamma = matrix(rnorm(L * K, 0, 0.1), nrow = L),
    D = abs(rnorm(K, 1, 0.1)) + 0.1,
    LKJ = diag(K) 
  )
}

# run the sampler
post_multilevel_nonreg_hmc <- mod_non_reg$sample(
  data = data_list,
  seed = 123,
  chains = 3,
  parallel_chains = 3,
  thin = 1,
  iter_warmup = 500,
  iter_sampling = 1000,
  init = init_fun,
  save_warmup = TRUE)

# save the posterior matrix
post_matrix_multilevel_nonreg_hmc <- post_multilevel_nonreg_hmc$draws(format="matrix")
save(post_matrix_multilevel_nonreg_hmc, file = "Task3/Posteriors/hmc_multilevel_noreg.RData")

# run the metropolis hastings algorithm for the non-regularized model

# load the model function
source("Task3/Stan_MH_Files/multilevel_logistic_mh.R")

# tune the standard deviation to get good results
proposal_sds <- c(0.01, 0.025, 0.05)

for (sd in proposal_sds) {
  result <- mh(X_multilevel_scaled, U_scaled, y_multilevel, country, n_samples = 1500,
               warmup = 500, proposal_sd = sd)
  
  cat("Chain with proposal sd:", sd, "\n")
  check_performance_mh(result)
  cat("----------------------------------------------\n")
}

# run three chains of the model with the tuned sd
post_multilevel_nonreg_mh <- bind_mh_multilevel(n_chains = 3, X = X_multilevel_scaled,
                                                U = U_scaled, y = y_multilevel,
                                                country, n_samples = 1500,
                                            warmup = 500, proposal_sd = 0.03)

# save the model
save(post_multilevel_nonreg_mh, file = "Task3/Posteriors/mh_multilevel_noreg.RData")


################### Run Regularized Multilevel Models ##########################

# compile stan model
mod_reg <- cmdstan_model("Task3/Stan_MH_Files/multilevel_logistic_regularized.stan")

# run the sampler
post_multilevel_reg_hmc <- mod_reg$sample(
  data = data_list,
  seed = 123,
  chains = 3,
  parallel_chains = 3,
  thin = 1,
  iter_warmup = 500,
  iter_sampling = 1000,
  init = init_fun,
  save_warmup = TRUE)

# save the posterior matrix
post_matrix_multilevel_reg_hmc <- post_multilevel_reg_hmc$draws(format="matrix")
save(post_matrix_multilevel_reg_hmc, file = "Task3/Posteriors/hmc_multilevel_reg.RData")

# run the metropolis hastings algorithm for the regularized model

# load the model function
source("Task3/Stan_MH_Files/multilevel_logistic_regularized_mh.R")

# tune the standard deviation to get good results
proposal_sds <- c(0.01, 0.025, 0.1)

for (sd in proposal_sds) {
  result <- mh(X_multilevel_scaled, U_scaled, y_multilevel, country, n_samples = 1500,
               warmup = 500, proposal_sd = sd)
  
  cat("Chain with proposal sd:", sd, "\n")
  check_performance_mh(result)
  cat("----------------------------------------------\n")
}

# run three chains of the model with the tuned sd
post_multilevel_reg_mh <- bind_mh_multilevel(n_chains = 3, X = X_multilevel_scaled,
                                                U = U_scaled, y = y_multilevel,
                                                country, n_samples = 1500,
                                                warmup = 500, proposal_sd = 0.025)

# save the model
save(post_multilevel_reg_mh, file = "Task3/Posteriors/mh_multilevel_reg.RData")


######################### Check Convergence ####################################

# from posteriors extract only log-likelihood
ll_hmc_noreg <- extract_matrix_hmc(post_multilevel_nonreg_hmc, 1, "ll")
ll_mh_noreg <- post_multilevel_nonreg_mh$ll
ll_hmc_reg <- extract_matrix_hmc(post_multilevel_reg_hmc, 1, "ll")
ll_mh_reg <- post_multilevel_reg_mh$ll

# arrange all log-likelihoods in a list
model_list <- list(ll_hmc_noreg, ll_mh_noreg, ll_hmc_reg, ll_mh_reg)

# set model names for plotting
model_names <- c("HMC - Weakly Informative", "MH - Weakly Informative",
                 "HMC - Regularized", "MH - Regularized")

# loop over all models
for (i in seq_along(model_list)) {
  
  # extract the ll as matrix and set respective column name
  ll_matrix <- matrix(model_list[[i]], ncol = 1)
  colnames(ll_matrix) <- "Log-Likelihood"
  
  # export each traceplot individually
  pdf(paste0("Task3/Figures/multilevel_traceplot_", model_names[[i]], ".pdf"),
      width = 6, height = 4)
  traceplots(ll_matrix, n_chains = 3, iterations = 1000)
  dev.off()
}


########################## Posterior Summaries #################################

# get the beta draws for all models
beta_draws_hmc_noreg <- post_matrix_multilevel_nonreg_hmc[
  ,grep("beta", colnames(post_matrix_multilevel_nonreg_hmc)), drop = FALSE]
beta_draws_hmc_reg <- post_matrix_multilevel_reg_hmc[
  ,grep("beta", colnames(post_matrix_multilevel_reg_hmc)), drop = FALSE]
beta_draws_mh_noreg <- post_multilevel_nonreg_mh$beta
beta_draws_mh_reg <- post_multilevel_reg_mh$beta

# define vector of individual-level predictor variables
var_names <- c("intercept", "satifaction_gov", "age", "unemployment")

# create long summaries of posterior matrix
beta_summary_hmc_noreg <- create_long_summary(beta_draws_hmc_noreg, "beta",
                                              countries, var_names)
beta_summary_hmc_reg <- create_long_summary(beta_draws_hmc_reg, "beta",
                                              countries, var_names)
beta_summary_mh_noreg <- create_long_summary_mh(beta_draws_mh_noreg, countries,
                                                var_names)
beta_summary_mh_reg <- create_long_summary_mh(beta_draws_mh_reg, countries,
                                              var_names)

# add information about the model type
beta_summary_hmc_noreg$model_type <- "Weakly Informative"
beta_summary_hmc_reg$model_type <- "Regularized"
beta_summary_mh_noreg$model_type <- "Weakly Informative"
beta_summary_mh_reg$model_type <- "Regularized"

# bind non-regularized and regularized models together
hmc_summary_beta <- rbind(beta_summary_hmc_noreg, beta_summary_hmc_reg)
mh_summary_beta <- rbind(beta_summary_mh_noreg, beta_summary_mh_reg)

# plot point estimates with confidence intervals and export the plot
# do this for both sampling algorithms
pdf("Task3/Figures/point_ci_multilevel_hmc.pdf", width = 10, height = 5)
plot_estimates_ci(hmc_summary_beta, countries, var_names)
dev.off()
pdf("Task3/Figures/point_ci_multilevel_mh.pdf", width = 10, height = 5)
plot_estimates_ci(mh_summary_beta, countries, var_names)
dev.off()


########################## Gamma Summary #######################################

# get the gamma draws for the hmc models
gamma_draws_hmc_noreg <- post_matrix_multilevel_nonreg_hmc[
  ,grep("gamma", colnames(post_matrix_multilevel_nonreg_hmc)), drop = FALSE]
gamma_draws_hmc_reg <- post_matrix_multilevel_reg_hmc[
  ,grep("gamma", colnames(post_matrix_multilevel_reg_hmc)), drop = FALSE]

# specify the names of group-level and individual parameters
group_var_names <- c("Intercept", "Corruption Control",
                     "Youth Unemployment", "Social Expenditures")
indiv_var_names <- c("Intercept", "Satisfaction Gov",
                     "Age", "Unemployment")

# compute gamma summary on the hmc posteriors
gamma_summary_hmc_noreg <- create_gamma_summary(gamma_draws_hmc_noreg,
                                                group_var_names, indiv_var_names)
gamma_summary_hmc_reg <- create_gamma_summary(gamma_draws_hmc_reg,
                                              group_var_names, indiv_var_names)

# export the tables as tex objects
export_table(gamma_summary_hmc_noreg, "Task3/Report/gamma_summary_hmc_noreg.tex")
export_table(gamma_summary_hmc_reg, "Task3/Report/gamma_summary_hmc_reg.tex")
