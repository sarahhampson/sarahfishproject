# Function to create a list of time series that qualify for novelty analysis
# Requirements are that there must be at least 10 unique, consecutive years in the dataset, and they must have greater than 2 species

select_timeseries <- function(data.c){
  # Want to make sure each time series is individual and unique
  #data.c = survey_data
  #bin_width = 1
  #length = 10 # years
  
  #### NON TIDYVERSE OPTION ####
  #Create frequency table of years for each timeseries
  #year_count <- data.table(data.c)[, .(number_of_years = length(unique(data.c$Year))), by = data.c$TimeSeriesID]
  #Remove timeseries with less than X no. years
  #year_count <- year_count[!(number_of_years<10)]
  
  #### TIDYVERSE OPTION ####
  year_count <- data.c %>% group_by(TimeSeriesID) %>%
    summarise(number_of_years = length(unique(Year)), first_year = min(Year), last_year = max(Year),
              first_to_last = last_year - first_year, full_data = (number_of_years == first_to_last)) %>%
  filter(full_data == TRUE)
  
  timeseries_list <- list(year_count$TimeSeriesID)
  
  return(timeseries_list)
}

