# Title: Basic results diagnostics/summaries
# Author: Sarah Hampson

# 0.  Load necessary dataframes -------------------------------------------

# Load full results csv
full_results_df <- read.csv("./detection results/full_results_df.csv")

# 1.  How many times did novelty occur? -----------------------------------

# Count of times novelty occurred
sum(full_results_df$taxnovel=="TRUE")
sum(full_results_df$funcnovel=="TRUE")

# Summary table for tax novelty
full_results_df %>%
  filter(taxnovel) %>%
  group_by(Site) %>%
  summarize(taxnovel_count = n()) %>%
  group_by(taxnovel_count) %>%
  summarize(site_count = n())
# Summary table for func novelty
full_results_df %>%
  filter(funcnovel) %>%
  group_by(Site) %>%
  summarize(taxnovel_count = n()) %>%
  group_by(taxnovel_count) %>%
  summarize(site_count = n())





# 2.  What country/biorealm did they occur in? ----------------------------

# Countries
# Summary table for tax novelty
full_results_df %>%
  filter(taxnovel) %>%
  group_by(Country) %>%
  summarize(taxnovel_count = n())
# Summary table for func novelty
full_results_df %>%
  filter(funcnovel) %>%
  group_by(Country) %>%
  summarize(funcnovel_count = n())

# BioRealm
# Summary table for tax novelty
full_results_df %>%
  filter(taxnovel) %>%
  group_by(BioRealm) %>%
  summarize(taxnovel_count = n())
# Summary table for func novelty
full_results_df %>%
  filter(funcnovel) %>%
  group_by(BioRealm) %>%
  summarize(funcnovel_count = n())
