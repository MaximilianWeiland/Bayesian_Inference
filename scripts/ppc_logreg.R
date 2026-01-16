# extract the replicated data as matrices
yrep_logreg_noreg_hmc <- post_matrix_logreg_nonreg_hmc[
  , grep("^y_rep", colnames(post_matrix_logreg_nonreg_hmc))]
yrep_logreg_noreg_mh <- post_logreg_nonreg_mh$y_rep
yrep_logreg_reg_hmc <- post_matrix_logreg_reg_hmc[
  , grep("^y_rep", colnames(post_matrix_logreg_reg_hmc))]
yrep_logreg_reg_mh <- post_logreg_reg_mh$y_rep

# wrap them inside a list
rep_matrices <- list(
  hmc_noreg = yrep_logreg_noreg_hmc,
  mh_noreg = yrep_logreg_noreg_mh,
  hmc_reg = yrep_logreg_reg_hmc,
  mh_reg = yrep_logreg_reg_mh
)

# get the actual turnout
y_turnout <- mean(y_germany)

# get replicate turnout and calculate p-values
p_values <- sapply(rep_matrices, function(mat) {
  rep_turnout <- rowMeans(mat)
  mean(rep_turnout >= y_turnout)
})

# store p-values inside a df and export
p_df_logreg <- data.frame(
  model = c("HMC - Weakly Informative", "MH - Weakly Informative",
            "HMC - Regularized", "MH - Regularized"),
  p_value = as.numeric(p_values)
)
export_table(p_df_logreg, "report/bayesian_p_values_logreg.tex")
