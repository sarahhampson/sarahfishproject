# Title: Novelty Detection Setup
# Author: Sarah Hampson


# 1.  Clear environment and install packages/functions --------------------

# Clear environment
rm(list = ls())
graphics.off()

# Set working directory
#setwd("YOUR OWN CHOICE HERE")
setwd("~/Documents/UNI/2023/Honours/BIOL6502/Data/sarahfish")

# Obtain required functions from 'functions' subfolder
sapply(list.files("./Functions", pattern="\\.R", full.names=TRUE), source)

#### Load packages
install.packages(c("mgcv", "vegan", "lme4", "nlme",
                   "DHARMa", "merTools", "shape",
                   "multcomp", "maptools", "sp", 
                   "divDyn", "plotrix", "raster",
                   "rgeos", "fun", "analogue",
                   "brms", "data.table", "dplyr",
                   "tidyverse", "funrar", "rfishbase",
                   "data.table", "svMisc", "data.table",
                   "tidyr"))

# Libraries required
library(tidyverse)
library(dplyr)

# 2.  Load time series data -----------------------------------------------

# RivFishTIME data
time_series_data <-read.csv("inputs/1873_2_RivFishTIME_TimeseriesTable.csv")
survey_data <- read.csv("inputs/1873_2_RivFishTIME_SurveyTable.csv")

# A little function to add dates to output files
date.wrap <- function(string, ext){
  paste0(string, " ", Sys.Date(), ext)
}

# 3.  Load fishmorph data -------------------------------------------------


# FishMORPH data plus change species variable name
fishmorph_data <- read.csv("inputs/FISHMORPH_Database.csv", sep = ";")
fishmorph_data <- fishmorph_data %>% dplyr::rename("Species" = "Genus.species")

# Align all species names across fishmorph with fishbase
fishmorph_data <- standardise.species.names(fishmorph_data)


# 4.  Select time series data ---------------------------------------------

# Function which gives filtered df for survey data and time series data
select.timeseries <- function(survey_data, time_series_data){
  
  # Remove extra spaces from the quarter column
  survey_data <- survey_data %>%
    mutate(Quarter = str_replace_all(Quarter, " ", ""))
  
  # If any 1/2 quarters appear replace the '/' with a '-'
  survey_data$Quarter <- gsub('/', '-', survey_data$Quarter)
  
  # Table of timeseries with min 10 consecutive yearly surveys and 2 species
  ts_table <- survey_data %>% group_by(TimeSeriesID) %>%
    dplyr::summarise(year_count = length(unique(Year)), first_year = min(Year), 
                     last_year = max(Year), first_to_last = last_year - first_year, 
                     full_data = (year_count == first_to_last+1), 
                     species_count = length(unique(Species)), 
                     quarter_count = length(unique(Quarter))) %>%
    filter(year_count >=10, species_count > 1) # removed filter full_data == TRUE
  
  # Modify survey data to only include timeseries from ts_table
  survey_data <- survey_data[survey_data$TimeSeriesID %in% 
                               ts_table$TimeSeriesID,]
  
  # Many TS are sampled in more than 1 quarter so I will filter out each quarter 
  # per time series (i.e. nested Quarter TS in one site)
  ts_q_table <- survey_data %>% group_by(TimeSeriesID, Quarter) %>%
    dplyr::summarise(year_count = length(unique(Year)),
                     species_count = length(unique(Species))) %>%
    filter(year_count >=10, species_count > 1)
  
  # Modify survey data to only include timeseries + quarters from ts_q_table
  survey_data <- merge(survey_data, ts_q_table, 
                       by.x = c("TimeSeriesID", "Quarter"))
  
  # Check quarter count frequency for the TS
  ts_data_qcount <- survey_data %>% group_by(TimeSeriesID) %>%
    dplyr::summarise(qcount = n_distinct(Quarter))
  
  # Check for frequency of TS with multiple surveys per quarter
  survey_quarter_table <- survey_data %>% group_by(TimeSeriesID, Year, Quarter) %>%
    dplyr::summarise(surveycount = length(unique(SurveyID)))
  # --> therefore I need to use the mean for abundance in quarterly time bins
  
  # Add site quarter index column
  survey_data <- survey_data %>% 
    mutate(sqID = paste0(TimeSeriesID, " Q", Quarter))
  
  return(survey_data)
}

# Select time series data using function
survey_data <- select.timeseries(survey_data)

# Modify time series to only include time series from new ts_table
ts_data <- time_series_data[time_series_data$TimeSeriesID %in%
                              survey_data$TimeSeriesID,]

# Standardise species names
survey_data <- standardise.species.names(survey_data)



# 5.  Make relative abundance dataframes ----------------------------------

# Create a list of relative abundance matrices using matrix maker function
# Note these are split into quarter time series
taxonomic_matrix_list <- seasonal.matrix.maker(survey_data, ts_data)


# Note --------------------------------------------------------------------

# Be sure to have taxonomic_matrix_list in your environment for both tax
# and func novelty detection

# 1.  Subset data for species in fishmorph --------------------------------

# Create a list of species from the selected TS
species_list <- survey_data %>% dplyr::summarise(Species = unique(Species))

# Create list of species missing from fishmorph
missing_species <- species_list %>%
  filter(!Species %in% fishmorph_data$Species)

# Identify TS which have these missing species
survey_data_missing_sp <- survey_data %>%
  filter(Species %in% missing_species$Species) %>% 
  dplyr::summarise(sqID = unique(sqID))

# Remove TS which do not have data from survey data
survey_data <- survey_data %>% 
  filter(!sqID %in% survey_data_missing_sp$sqID)

# Remove sites which do not appear in the new survey data
ts_data <- ts_data %>% 
  filter(TimeSeriesID %in% survey_data$TimeSeriesID)

# Make new species list
species_list <- survey_data %>% dplyr::summarise(Species = unique(Species))

# Cleanup
rm(missing_species, survey_data_missing_sp)