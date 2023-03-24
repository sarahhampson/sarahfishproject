# Title: Trait data setup
# Author: Sarah Hampson


# 1.  Identify species ----------------------------------------------------

# Import survey_data if not already in environment
survey_data<- read.csv("./inputs/survey_data.csv")

# Make list of unique species
species_list <- survey_data %>% dplyr::summarise(Species = unique(Species))

# 2.  Fishmorph traits ----------------------------------------------------

# Make a list of traits that we are interested in from fishmorph
traitlist_FM <- as.list(c("MBl", "BEl", "VEp", "REs", "OGp", "RMl", "BLs", "PFv", "PFs", "CPt"))

# Use fishmorph_trait_finder to select all data from fishmorph
trait_data_FM <- fishmorph_trait_finder(traitlist_FM, fishmorph_data, species_list)

# 3.  Trophic level -------------------------------------------------------

# Retrieve trophic level data from rfishbase
trophic_data_FB <- ecology(species_list_v2$Species, 
                           fields = c("Species", "FoodTroph", "DietTroph"))

# Make list of species missing from troph data
missing_troph_sp <- species_list_v2 %>% 
  filter(!Species %in% trophic_data_FB$Species)

# Show more missing species from above troph list and filter from df
missing_troph_sp2 <- trophic_data_FB %>% filter(is.na(FoodTroph), is.na(DietTroph))
trophic_data_FB <- trophic_data_FB %>% filter(!Species %in% missing_troph_sp2$Species)

# Make sure no double ups in troph data from FB by taking mean diet and food troph
trophic_data_FB <- trophic_data_FB %>% 
  group_by(Species) %>% 
  dplyr::summarise(FoodTroph = mean(FoodTroph, na.rm=T),
                   DietTroph = mean(DietTroph, na.rm=T))

# There are 151+37 so 188 species missing, which I found information on via online fishbase

# Read in csv of model based fishbase species without specific trophic information
troph_Est_data <- read_csv('./inputs/trophic_level_FBModels.csv', col_names = T)
troph_Est_data$Species <- gsub("_", " ", troph_Est_data$Species)
troph_Est_data <- troph_Est_data %>% 
  rename(troph_level = Trophic_level) %>% 
  mutate(troph_cat = "Est") %>% 
  dplyr::select(-Cat)

# Filter troph diet data first
troph_data_diet <- trophic_data_FB %>% 
  filter(!is.na(DietTroph)) %>% 
  rename(troph_level = DietTroph) %>% 
  dplyr::select(-FoodTroph) %>% 
  mutate(troph_cat = "Diet")

# Filter troph food data
troph_data_food <- trophic_data_FB %>% 
  filter(!Species %in% troph_data_diet$Species) %>% 
  rename(troph_level = FoodTroph) %>% 
  dplyr::select(-DietTroph) %>% 
  mutate(troph_cat = ifelse(troph_level>0, "Food", NA))

# Combine the data into one trophic df
troph_data <- rbind(troph_Est_data, troph_data_diet, troph_data_food)

# Cleanup
rm(troph_data_diet, troph_data_food, troph_Est_data, trophic_data_FB,
   missing_troph_sp, missing_troph_sp2)

# 4.  Testing for collinearity --------------------------------------------

# Checking to see whether co-linearity exists between traits
# Using VIF (variance inflation factor)
# If VIF = 1 there is no correlation
# If 1<VIF<5 there is moderate correlation
# If VIF> 5 there is strong correlation

# Install VIF package and define multiple linear regression model
install.packages("car")
library(car)
# Trial all versions of this where each trait is the response
colin.model <- lm(CPt ~ MBl + BEl + VEp + REs + OGp + RMl + BLs + PFv + PFs + troph_level, data=full_trait_data)
vif(colin.model)
summary(colin.model)
plot(colin.model)

# Cleanup
rm(colin.model)


# 6.  Make functional dissimilarity matrix --------------------------------

# Need to make a multitrait dissimilarity matrix between all species using gowdis
func_dissims_matrix <- as.matrix(gowdis(full_trait_data[,2:12]), dimnames = list(full_trait_data$Species))
dimnames(func_dissims_matrix) <- list(full_trait_data$Species, full_trait_data$Species)
attributes(func_dissims_matrix) # All traits should be numeric

# Change to numeric
class(func_dissims_matrix) <- "numeric"

# See first 5 rows and cols of diss matrix
func_dissims_matrix[504:509, 504:509]

# 5.  Save trait and dissims data -----------------------------------------------------

# Put all trait data together and remove unnecessary columns
full_trait_data <- left_join(trait_data_FM, troph_data, by="Species")
full_trait_data <- full_trait_data %>% dplyr::select(-troph_cat)
full_trait_data$troph_level <- as.numeric(full_trait_data$troph_level)

# Create filename and write csv to inputs for later use
filename_traits <- date.wrap("./inputs/trait_data", ".csv")
write.csv(full_trait_data, filename_traits)

# Save func dissims matrix as csv to inputs for later use
filename_dissims <- date.wrap("./inputs/func_dissims_matrix", ".csv")
write.csv(func_dissims_matrix, filename_dissims)

# Cleanup
rm(filename_dissims, filename_traits)
