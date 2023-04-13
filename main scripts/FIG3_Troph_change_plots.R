# Title: Plots for trophic level changes
# Author: Sarah Hampson


# 0.  Load necessary dataframes -------------------------------------------

# All community size change data
comms_troph_df <- read.csv("./inputs/model_data/trait_models/troph_df.csv")

# Taxonomic troph change model summary df
tax_troph_plot_df <- read.csv("./outputs/troph/tax_lm_troph.csv")
tax_Abtroph_plot_df <- read.csv("./outputs/troph/tax_lm_Abtroph.csv")
tax_sdtroph_plot_df <- read.csv("./outputs/troph/tax_lm_sdtroph.csv")
func_troph_plot_df <- read.csv("./outputs/troph/func_lm_troph.csv")
func_Abtroph_plot_df <- read.csv("./outputs/troph/func_lm_Abtroph.csv")
func_sdtroph_plot_df <- read.csv("./outputs/troph/func_lm_sdtroph.csv")

# 1.  Subset and clean dataframes ---------------------------------------------------

# Take only intercept and novel=True estimates
tax_troph_plot_df <- tax_troph_plot_df[1:2,]
tax_Abtroph_plot_df <- tax_Abtroph_plot_df[1:2,]
tax_sdtroph_plot_df <- tax_sdtroph_plot_df[1:2,]
func_troph_plot_df <- func_troph_plot_df[1:2,]
func_Abtroph_plot_df <- func_Abtroph_plot_df[1:2,]
func_sdtroph_plot_df <- func_sdtroph_plot_df[1:2,]

# Absolute change dataframes need to be exponentially transformed
# but only transform the mean estimate and intervals
tax_Abtroph_plot_df <- exp(tax_Abtroph_plot_df[, 3:7])
func_Abtroph_plot_df <- exp(func_Abtroph_plot_df[, 3:7])

# Put in a new "cat" variable for "Non novel" and "Novel" rows
tax_troph_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_Abtroph_plot_df$taxnovel <- c("FALSE", "TRUE")
tax_sdtroph_plot_df$taxnovel <- c("FALSE", "TRUE")
func_troph_plot_df$funcnovel <- c("FALSE", "TRUE")
func_Abtroph_plot_df$funcnovel <- c("FALSE", "TRUE")
func_sdtroph_plot_df$funcnovel <- c("FALSE", "TRUE")

# 2.  Absolute troph change plots ------------------------------------------

# Taxonomic
ggplot(data = tax_Abtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=taxnovel, Abtraitchange), trim=T) +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=taxnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Taxonomic novelty")) +
  labs(x="Community category", y=expression("|" *Delta*" community mean troph |"))
# Save plot
ggsave(
  filename="taxAbtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/troph",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# Functional
ggplot(data = func_Abtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=funcnovel, Abtraitchange), trim=T) +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=funcnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Functional novelty")) +
  labs(x="Community category", y=expression("|" *Delta*" community mean troph |"))
# Save plot
ggsave(
  filename="funcAbtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# 3.  Raw troph change plots -----------------------------------------------

# Taxonomic
ggplot(data = tax_troph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=taxnovel, traitchange), trim=F) +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=taxnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Taxonomic novelty")) +
  labs(x="Community category", y=expression(Delta*" community mean troph"))
# Save plot
ggsave(
  filename="taxrawtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# Functional
ggplot(data = func_troph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=funcnovel, traitchange), trim=F) +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=funcnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Functional novelty")) +
  labs(x="Community category", y=expression(Delta*" community mean troph"))
# Save plot
ggsave(
  filename="funcrawtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# 4.  Standard deviation change plots -------------------------------------

# Taxonomic
ggplot(data = tax_sdtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=taxnovel, sdchange), trim=F) +
  geom_point(aes(taxnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(taxnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=taxnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Taxonomic novelty")) +
  labs(x="Community category", y=expression(Delta* " standard deviation in troph"))
# Save plot
ggsave(
  filename="taxsdtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# Functional
ggplot(data = func_sdtroph_plot_df) + 
  geom_violin(data=comms_troph_df, aes(x=funcnovel, sdchange), trim=F) +
  geom_point(aes(funcnovel, y=Mean_Estimate)) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=CI_Lower, ymax=CI_Upper), width=0.15, linetype=1) +
  geom_errorbar(aes(as.factor(funcnovel), ymin=PI_Lower, ymax=PI_Upper), width=0.15, linetype=2) +
  #geom_jitter(data=comm_troph_df, aes(x=funcnovel, y=traitchange), width=0.2, alpha=0.25) 
  theme_classic() +
  theme(axis.ticks.x=element_blank()) +
  scale_x_discrete(labels = c("Not novel", "Functional novelty")) +
  labs(x="Community category", y=expression(Delta*" standard deviation in troph"))
# Save plot
ggsave(
  filename="funcsdtrophplot", plot = last_plot(),
  device = "pdf", path = "./figures/MBl",
  scale = 1, width = 10, height = 15,
  units = "cm", dpi = 300)

# 5.  Cleanup -------------------------------------------------------------

rm(comms_troph_df, func_Abtroph_plot_df, func_glm_site_df, func_sdtroph_plot_df,
   func_troph_plot_df, tax_Abtroph_plot_df, tax_glm_site_df, tax_sdtroph_plot_df,
   tax_troph_plot_df)
