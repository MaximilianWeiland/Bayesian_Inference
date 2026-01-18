# Drivers of Political Participation - Bayesian Analysis
This project investigates the research question of what factors motivate people to participate in national elections. Hypotheses are drawn for personal motivations as well as for the interplay between these individual-level drivers and structural country-level factors. The expectations are tested with Bayesian logistic and multilevel logistic regressions. The Bayesian MCMC sampling is achieved via model formulation in Stan files as well as a manual implementation with the Metropolis Hastings algorithm.

## Quick Run
To run the project, simply execute the `master` script. This script loads custom functions from `functions` and executes all individual R scripts collected under [scripts](scripts).

## Structure

### data
Data from the European Social Survey (ESS) is used to measure the degree of political participation. Information on country-level degrees of social expenditures, corruption levels and youth unemployment is obtained by the OECD data portal and the World Bank.

### figures
All plots created for the report saved as png-files.

### model files
All models are implemented as Stan files which are stored in this folder. Additionally, to get a better understanding of the sampling process, each model is also implemented via the Metropolis Hastings algorithm in an additional R script.

### posteriors
Contains the posterior matrices for all models.

### report
Here you can find the final analysis report as a PDF-file.

### scripts
All R scripts for data preprocessing, model loading, diagnostics and visualizations. These files get loaded when running the `master` file.
