# load all created functions
source("functions.R")

# call the data preparation script
source("scripts/data_prep.R")

# run all logistic regression models, check convergence and create summary
source("scripts/exec_diagnostics_logreg.R")

# visualize marginal effects of all variables (population-average prediction)
source("scripts/prediction_visualization_logreg.R")

# conduct ppc of all logistic regression models
source("scripts/ppc_logreg.R")

# run all multilevel models, check convergence and create posterior summaries
source("scripts/exec_diagnostics_multilevel.R")

# conduct ppc of all multilevel models
source("scripts/ppc_multilevel.R")
