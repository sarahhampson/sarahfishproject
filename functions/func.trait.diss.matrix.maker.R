# A function which takes a trait dataframe and converts this into
# dissimilarity matrices

# Define function
func.diss.matrix.maker <- function(traitdata, surveydata, tsdata){
  
  # Verification variables
  # traitdf <- pilot_body_length_data
  
  # Convert trait dataframe to dissimilarity matrix
  trait.diss.mat <- lapply(traitdf[, -traitdf$Species],
    matrix(pilot_body_length_data, 
         dimnames = list(pilot_body_length_data$Species, "Length")))
  
  # Convert matrix to numeric form <- may need changing 
  # depending on traits of interest
  class(traitdf) <- mat.type
  
  
  matrix(pilot_body_length_data$Length, 
         dimnames = list(pilot_body_length_data$Species, "Length"))
  
  
  
}
