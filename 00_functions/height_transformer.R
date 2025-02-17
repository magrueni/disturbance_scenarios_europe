# =============================================================================
# Title: Height Prediction Model for Tree Species
# Author: Marc Grünig
# Date: 2022
# 
# Description:
# This script implements a height prediction model for various tree species 
# using yield table data and simulation data. The primary focus is to estimate 
# dominant heights (Oberhöhe) based on species-specific models, leveraging 
# relationships between maximum, mean, and minimum heights. The analysis includes:
# 
# 1. Loading and preprocessing yield table data (Ertragstafel) with German 
#    species names translated into species codes.
# 2. Visual and statistical exploration of height data, including error checks 
#    and visual comparisons.
# 3. Fitting and evaluating linear and mixed-effects models to predict dominant 
#    heights for different species.
# 4. Predicting heights for real-world simulation datasets from a SQLite database 
#    under various scenarios, including cases with missing mean or maximum heights.
# 5. Incorporating approaches to handle species combinations and dominant species 
#    scenarios using adjustment factors.
# 
# Dependencies:
# - R Packages: ggplot2, dplyr, lme4, MuMIn, DBI, RSQLite
# - Custom Functions: yield_table.R (for loading yield table data)
# - External Data: SQLite database (simulation_db_users_v4.sqlite)
# 
# Outputs:
# - Visualizations: Scatter plots and diagnostic plots of model predictions.
# - Model coefficients and predictions saved for further analysis.
# 
# Notes:
# - Ensure that all required packages are installed before running the script.
# - Update file paths for the yield table script, SQLite database, and output files 
#   as necessary.
# =============================================================================


# Load libraries
library(ggplot2)
library(dplyr)
library(lme4)
library(MuMIn)
library(DBI)
library(RSQLite)
library(Matrix)

# Load yield table
source("./publi/functions/yield_table.R")
etafel <- read.delim(text = etafel.str)

# Map species names to codes
etafel <- etafel %>%
  mutate(
    species = case_when(
      ET_BAUMART_NAME %in% c("Birke", "Esche, Ungarn", "europäische Pappel", "Hainbuche", "Hybridpappel", "Robinie, Ungarn", "Schwarzerle") ~ "Laub",
      ET_BAUMART_NAME %in% c("Kiefer", "Schwarzkiefer") ~ "Kiefer",
      TRUE ~ ET_BAUMART_NAME
    ),
    species_short = factor(
      species,
      labels = c("piab", "abal", "lade", "pisy", "fasy", "quro", "laub")
    )
  )

# Visual checks
ggplot(etafel, aes(x = ET_TAB_OBERHOEHE, y = H_max, col = species_short)) +
  geom_point(pch = ".") +
  geom_abline()


# Fit linear and mixed-effects models
hcorrmodel1 <- lm(ET_TAB_OBERHOEHE ~ hdom_est, data = etafel)
hcorrmodel_re <- lmer(ET_TAB_OBERHOEHE ~ hdom_est + (hdom_est | species_short), data = etafel)

# Save coefficients
coefs <- coef(hcorrmodel_re)$species_short
coefs_glob <- fixef(hcorrmodel_re)
all_coefs <- rbind(global = coefs_glob, coefs) %>% cbind(., rownames(.))
colnames(all_coefs) <- c("intercept", "slope", "species")

# Plot model predictions
etafel$oberh_pred <- predict(hcorrmodel_re, etafel)

ggplot(etafel, aes(x = ET_TAB_OBERHOEHE, y = oberh_pred, col = species_short)) +
  geom_point(pch = ".") +
  geom_abline()

# Load simulation data
sim_db <- dbConnect(SQLite(), "/.../forest_simulation_db_v1.sqlite")
sim_tables <- dbListTables(sim_db) %>% grep("simulationdata", ., value = TRUE)

simdat_all <- sim_tables %>%
  lapply(function(tab) {
    simdat <- dbReadTable(sim_db, tab)
    if (nrow(simdat) > 20000) {
      simdat <- simdat %>% sample_n(20000)
    }
    simdat %>% select(
      uniqueID, simulationID, Year, Species1, Proportion1, Species2, Proportion2,
      Species3, Proportion3, Species4, Proportion4, Species5, Proportion5,
      MeanHeight, MaxHeight, MinHeight
    )
  }) %>%
  bind_rows()

dbDisconnect(sim_db)

# Model-based predictions
newdat <- simdat_all %>%
  filter(!is.na(MaxHeight) & !is.na(MeanHeight)) %>%
  mutate(
    hdom_est = (MaxHeight + MeanHeight) / 2,
    hdom_est2 = (((MeanHeight * 2) - MinHeight) + MeanHeight) / 2,
    species_short = case_when(
      Species1 %in% c("bepe", "frex", "potr", "cabe", "poca", "rops", "algl") ~ "laub",
      !Species1 %in% c("piab", "abal", "lade", "pisy", "fasy", "quro", "laub") ~ "other",
      TRUE ~ Species1
    ),
    pred_oberh = predict(hcorrmodel_re, ., allow.new.levels = TRUE)
  ) %>%
  mutate(
    pred_oberh = ifelse(pred_oberh / MaxHeight < 0.8, MaxHeight * 0.8, pred_oberh),
    pred_oberh = ifelse(pred_oberh > MaxHeight, MaxHeight, pred_oberh)
  )

# Plot predictions
ggplot(newdat, aes(x = MaxHeight, y = pred_oberh, col = species_short)) +
  geom_point(pch = ".") +
  geom_abline()

# Handle cases with missing MeanHeight
sim_nw <- simdat_all %>%
  filter(is.na(MeanHeight) & !is.na(MaxHeight)) %>%
  mutate(
    MeanHeight = (MinHeight + MaxHeight) / 2,
    hdom_est = (MaxHeight + MeanHeight) / 2
  )

# Species factor-based adjustment for predictions
newdat <- sim_nw %>%
  mutate(
    species_short = case_when(
      Species1 %in% c("bepe", "frex", "potr", "cabe", "poca", "rops", "algl") ~ "laub",
      !Species1 %in% c("piab", "abal", "lade", "pisy", "fasy", "quro", "laub") ~ "other",
      TRUE ~ Species1
    ),
    pred_oberh_method2 = predict(hcorrmodel_re, ., allow.new.levels = TRUE)
  )

# Species-specific scaling factors for height
newdat <- sim_nw %>%
  mutate(
    pred_oberh_method4 = case_when(
      Species1 %in% c("piab", "abal", "lade") ~ MaxHeight / 1.138,
      Species1 %in% c("pisy") ~ MaxHeight / 1.18,
      Species1 %in% c("fasy") ~ MaxHeight / 1.132,
      Species1 %in% c("quro") ~ MaxHeight / 1.184,
      TRUE ~ MaxHeight / 1.14
    )
  )

# Final plot
ggplot(newdat, aes(x = MaxHeight, y = pred_oberh_method4, col = Species1)) +
  geom_point(pch = ".", size = 3) +
  geom_abline()
