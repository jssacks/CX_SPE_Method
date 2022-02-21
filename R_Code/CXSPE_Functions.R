####
####
####
library(tidyverse)
library(readr)


###Function for importing skyline output csv files, removing unneeded information
# and renaming columns
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