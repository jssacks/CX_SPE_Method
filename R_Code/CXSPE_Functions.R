####
####
####
library(tidyverse)
library(readr)


###Function for importing skyline output csv files, removing unneeded information
# and renaming columns for analytical validation experiment
csv_parse <- function(d) {
  d1 <- read_csv(d, col_types = cols(Area = col_number()))
  d2 <- d1 %>%
    select('Replicate Name', 'Precursor Ion Name', 'Area') %>%
    rename(
      Rep = 'Replicate Name',
      Compound = 'Precursor Ion Name'
    )
  return(d2)
}

###Function for importing skyline output csv files, removing unneeded information
# and renaming columns for environmental samples experiment
ES_csv_parse <- function(d) {
  d1 <- read_csv(d, col_types = cols(Area = col_number()))
  d2 <- d1 %>%
    select('Replicate Name', 'Precursor Ion Name', 'Area', "Background", "Mass Error PPM") %>%
    rename(
      Rep = 'Replicate Name',
      Compound = 'Precursor Ion Name',
      Mass.Error.PPM = `Mass Error PPM`
    ) %>%
    mutate('Background' = as.numeric(Background),
           'Mass.Error.PPM' = as.numeric(Mass.Error.PPM))
  return(d2)
}
