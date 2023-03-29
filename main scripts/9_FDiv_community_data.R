# Title: PCA axes for morphological traits
# Author: Sarah Hampson

# 0.  Load necessary dataframes -------------------------------------------

# Load full results csv
trait_data <- read.csv("./inputs/data/trait_data.csv")

# Install mFD package
install.packages("mFD")
library(mFD)

# Get one species x abundance matrix
relA_data <- data.frame(taxonomic_matrix_list[[1]], check.names = FALSE)
relA_species <- colnames(relA_data)

# Drop the species from the trait_data that aren't in the abundance data
trait_data <- trait_data %>% filter(Species %in% relA_species)
rownames(trait_data) <- trait_data$Species
trait_data <- trait_data %>% dplyr::select(-X, -Species)


test <- dbFD(trait_data,
             relA_data,
             w.abun=T,
             calc.FRic=T,
             m=5,
             calc.CWM=T,
             calc.FDiv = T,
             print.pco=T,
             )

traits_cat <-data.frame(matrix(ncol=2, nrow=11))
colnames(traits_cat) <- c("trait_name", "trait_type")
traits_cat$trait_name <- names(trait_data)
traits_cat$trait_type <- "Q"


# Make distrance matrix
testmFD <- mFD::funct.dist(
  sp_tr         = trait_data,
  tr_cat        = traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)
