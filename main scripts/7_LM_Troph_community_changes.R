# Title: LM models for Trophic level changes in community mean and sd
# Author: Sarah Hampson

# Note --> trophic level is a continuous and normally distributed trait

# 0.  Load necessary dataframes -------------------------------------------

# Read results csvs
#tax_results_df <- read.csv("./detection results/tax_results_df.csv")
#func_results_df <- read.csv("./detection results/func_results_df.csv")
full_results_df <- read.csv("./detection results/full_results_df.csv")
trait_data <- read.csv("./inputs/data/trait_data.csv")

# Need taxonomic matrix list, if not already in environment run this
survey_data <- read.csv("./inputs/data/survey_data.csv")
ts_data <- read.csv("./inputs/data/ts_data.csv")
taxonomic_matrix_list <- seasonal.matrix.maker(survey_data, ts_data)

# 1.  Retrieve community trophic level data  ---------------------------------------

# Includes:
# absolute change in community mean from one time point to the next
# raw change in community mean from one time point to the next
# absolute change in standard deviation from one time point to the next
# raw change in standard deviation from one time point to the next

# All of these variables are abundance weighted!

# Make list of communities and include year it emerges and the year of data beforehand
comms <- full_results_df %>% 
  mutate(bins.before = as.numeric(bins)-as.numeric(bin.lag)) %>% 
  dplyr::select(sqID, bins, bins.before)

# Abundance weighted community mean trait
comms_troph_data_list <- find.commS.trait.data(comms, taxonomic_matrix_list, trait_data, "all")

# 2.  Merge with detection results ----------------------------------------

# Row bind these into tax and func comm size data
comms_troph_data <- do.call("rbind", comms_troph_data_list)
comms_troph_df <- data.frame(comms_troph_data)

# Change character columns to numeric
comms_troph_df[, 2:9] <- lapply(2:9, function(x){
  as.numeric(comms_troph_df[[x]])
})

# Merge with full results df
full_results_df_troph <- merge(full_results_df, comms_troph_df,
                              by.y=c("sqID", "YearofComm"),
                              by.x=c("sqID", "bins"))

# 3.  Save community troph dataframe ---------------------------------------

# Save community troph data into the model_data -> trait_models folder
write.csv(full_results_df_troph, "./inputs/model_data/trait_models/troph_df.csv")

# 4.  Add and scale fixed effects -----------------------------------------

temp_df <- lapply(split(full_results_df_troph, f=full_results_df_troph$sqID),
                  function(x){
                    x$TS_Length <- nrow(x)
                    x$Yrs_passed <- x$bins - min(x$bins)
                    x$Yr_pos <- 1:nrow(x)
                    return(x)
                  })
full_results_df_troph <- do.call("rbind", temp_df)

# Scale these new variables
full_results_df_troph$bin.lagS <- scale(full_results_df_troph$bin.lag)[,1]
full_results_df_troph$Yrs_passedS <- scale(full_results_df_troph$Yrs_passed)[,1]
full_results_df_troph$logYr_posS <- scale(log(full_results_df_troph$Yr_pos))[,1]

# 5.  Troph changes ~ tax novelty models -----------------------------------

# Model trait change as a function of whether community is taxnovel or not
# Including same fixed/randomeffects as in other time point models
tax.lm.troph <- lm(traitchange ~ taxnovel + bin.lagS + logYr_posS, 
                  data = full_results_df_troph)

# I have to remove data which have absolute zero change
full_results_df_Abtroph <- full_results_df_troph %>% 
  filter(Abtraitchange>0)
# 0.0567 --> ~ 6% of datapoints had zero change in troph
# Model absolute trait change as a function of whether community is
# tax novel or not
tax.lm.Abtroph <- lm(log(Abtraitchange) ~ taxnovel + bin.lagS + logYr_posS, 
                    data = full_results_df_Abtroph)

# Model change in standard deviation of troph as a function of whether
# community is tax novel or not
tax.lm.sdtroph <- lm(sdchange ~ taxnovel + bin.lagS + logYr_posS, 
                    data = full_results_df_troph)
# See summary
summary(tax.lm.troph)
summary(tax.lm.Abtroph)
summary(tax.lm.sdtroph)

# 6.  Troph changes ~ func novelty models ----------------------------------

# Model trait change as a function of whether community is funcnovel or not
func.lm.troph <- lm(traitchange ~ funcnovel + bin.lagS + logYr_posS, 
                   data = full_results_df_troph)
func.lm.Abtroph <- lm(log(Abtraitchange) ~ funcnovel + bin.lagS + logYr_posS, 
                     data = full_results_df_Abtroph)
func.lm.sdtroph <- lm(sdchange ~ funcnovel + bin.lagS + logYr_posS, 
                     data = full_results_df_troph)
# See summary
summary(func.lm.troph)
summary(func.lm.Abtroph)
summary(func.lm.sdtroph)

# 7.  Create and save plot dataframes -------------------------------------

# Use function to create plot dataframes for each model
tax_troph_plot_df <- make.lm.plot.df(tax.lm.troph)
tax_Abtroph_plot_df <- make.lm.plot.df(tax.lm.Abtroph, log=T)
tax_sdtroph_plot_df <- make.lm.plot.df(tax.lm.sdtroph)
func_troph_plot_df <- make.lm.plot.df(func.lm.troph)
func_Abtroph_plot_df <- make.lm.plot.df(func.lm.Abtroph, log=T)
func_sdtroph_plot_df <- make.lm.plot.df(func.lm.sdtroph)

# Save the dataframes
write.csv(tax_troph_plot_df, "./outputs/troph/tax_lm_troph.csv")
write.csv(tax_Abtroph_plot_df, "./outputs/troph/tax_lm_Abtroph.csv")
write.csv(tax_sdtroph_plot_df, "./outputs/troph/tax_lm_sdtroph.csv")
write.csv(func_troph_plot_df, "./outputs/troph/func_lm_troph.csv")
write.csv(func_Abtroph_plot_df, "./outputs/troph/func_lm_Abtroph.csv")
write.csv(func_sdtroph_plot_df, "./outputs/troph/func_lm_sdtroph.csv")

# 8.  Cleanup ------------------------------------------------------------

rm(full_results_df, trait_data, taxonomic_matrix_list, survey_data, ts_data, 
   comms_troph_data_list, comms_troph_data, comms_troph_df, tax.lm.troph,
   tax.lm.Abtroph, tax.lm.sdtroph, tax_troph_plot_df, tax_Abtroph_plot_df,
   tax_sdtroph_plot_df, func.lm.troph, func.lm.Abtroph, func.lm.sdtroph,
   func_troph_plot_df, func_Abtroph_plot_df, func_sdtroph_plot_df,
   full_results_df_troph, full_results_df_Abtroph, comms, temp_df)








