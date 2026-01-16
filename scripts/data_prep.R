# load all required libraries
library(bayesplot)
library(xtable)
library(cmdstanr)
library(MASS)
library(abind)

# load the data
ess <- read.csv("data/ESS11.csv")
control_of_corruption <- read.csv("data/wb_control_of_corruption.csv")
youth_unemployment <- read.csv("data/wb_youth_unemployment.csv")
social_expenditures <- read.csv("data/oecd_social_expenditures.csv")

# remove unnecessary rows and columns of dataframes
# note: I needed to remove some more countries as they are not present in ess
control_of_corruption <- control_of_corruption[1:8, c(2, 5)]
control_of_corruption <- control_of_corruption[c(-5, -6), ]
youth_unemployment <- youth_unemployment[1:8, c(4, 5)]
youth_unemployment <- youth_unemployment[c(-7, -8), ]
social_expenditures <- social_expenditures[-7, c("REF_AREA", "OBS_VALUE")]
social_expenditures <- social_expenditures[-c(2, 8), ]

# define which countries to analyze (in multilevel analysis)
countries <- c("DE", "FR", "HU", "IT", "NO", "PL")

# define country remapping
country_map <- c(
  "DEU" = "DE",
  "FRA" = "FR",
  "HUN" = "HU",
  "ITA" = "IT",
  "POL" = "PL",
  "NOR" = "NO"
)

# apply remapping to all country-level dataframes
control_of_corruption$country_name <- country_map[control_of_corruption$Country.Code]
youth_unemployment$country_name <- country_map[youth_unemployment$Country.Code]
social_expenditures$country_name <- country_map[social_expenditures$REF_AREA]

# rename columns of covariate data
names(control_of_corruption)[2] <- "corruption_control"
names(youth_unemployment)[2] <- "share_youth_unemployment"
names(social_expenditures)[2] <- "social_exp_gdp"

# merge them to one country-level df
country_level_df <- merge(control_of_corruption, youth_unemployment,
                          by = "country_name", all.x = T)
country_level_df <- merge(country_level_df, social_expenditures,
                          by = "country_name", all.x = T)
country_level_df <- country_level_df[, c("country_name", "corruption_control",
                                         "share_youth_unemployment",
                                         "social_exp_gdp")]

# filter ess data for all countries to consider in multilevel analysis
ess_multilevel <- ess[ess$cntry %in% countries, ]

# join country data
ess_multilevel <- merge(ess_multilevel, country_level_df, by.x = "cntry",
                        by.y = "country_name", all.x = T)

ess_germany <- ess[ess$cntry == "DE", ]

# subset the df for observations that have valid entries in all columns
valid_germany <- with(ess_germany,
                      vote %in% c(1, 2) &
                        stfgov %in% 1:10 &
                        agea != 999 &
                        hincsrca %in% 1:8) 
ess_germany_subset <- ess_germany[valid_germany, ]

valid_multilevel <- with(ess_multilevel,
                      vote %in% c(1, 2) &
                        stfgov %in% 1:10 &
                        agea != 999 &
                        hincsrca %in% 1:8)
ess_multilevel_subset <- ess_multilevel[valid_multilevel, ]

# extract dependent variable for both models and convert to dummy (1 = yes)
y_germany <- ifelse(ess_germany_subset$vote == 1, 1, 0)
y_multilevel <- ifelse(ess_multilevel_subset$vote == 1, 1, 0)

# extract the grouping variable
country <- as.numeric(factor(ess_multilevel_subset$cntry))

# extract predictor variables at the individual level for both models
satisfaction_gov_germany <- ess_germany_subset$stfgov
satisfaction_gov_multilevel <- ess_multilevel_subset$stfgov
age_germany <- ess_germany_subset$agea
age_multilevel <- ess_multilevel_subset$agea
unemployment_germany <- ifelse(ess_germany_subset$hincsrca == 5, 1, 0)
unemployment_multilevel <- ifelse(ess_multilevel_subset$hincsrca == 5, 1, 0)

# extract country-level predictors for the multilevel model
corruption_control <- tapply(ess_multilevel_subset$corruption_control,
                             as.factor(country), min)
share_youth_unemployment <- tapply(ess_multilevel_subset$share_youth_unemployment,
                                   as.factor(country), min)
social_exp_gdp <- tapply(ess_multilevel_subset$social_exp_gdp,
                         as.factor(country), min)

# create data objects for both the individual and multilevel model
X_ger_scaled <- cbind(1, scale(cbind(satisfaction_gov_germany, age_germany)),
                      unemployment_germany)
X_multilevel_scaled <- cbind(1, scale(cbind(satisfaction_gov_multilevel,
                                            age_multilevel)),
                             unemployment_multilevel)
U_scaled <- cbind(1, corruption_control, scale(cbind(share_youth_unemployment,
                                 social_exp_gdp)))

# set dimensionality variables for the stan models
K <- ncol(X_multilevel_scaled)
N <- length(y_multilevel)
J <- max(country)
L <- ncol(U_scaled)
