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


#############Function for determininng R2 values from calibration curves in ALOHA seawater
Calibcurv_calc <- function(dataset, compound) {

####select for just a single compound
  dat <- dataset %>%
    filter(MF == compound)

#####Calculate average no spike value
  NS.av <- dat %>%
    filter(replicate == "NS") %>%
    summarise("NS.av" = mean(Adjusted_Area))

######Subtract out no-spike average value
  Dat.lin <- dat %>%
    filter(!replicate == "NS") %>%
    filter(!replicate == "SAf") %>%
    mutate("Adjusted_Area" = Adjusted_Area - NS.av$NS.av) %>%
    mutate(Spike.Conc = as.numeric(Spike.Conc))
  
#####Run linear model on data to determine R2 value
  ols.linmod <- lm(formula = Adjusted_Area ~ Spike.Conc,
                   data = Dat.lin)
  ols.r2 <- summary(ols.linmod)$r.squared
  
#####add R2 value to data frame
  dat.lin.output <- Dat.lin %>%
    mutate("R2" = ols.r2) %>%
    select(MF, column, Fraction, z, R2)

#remove duplicate value
  dat.lin.output <- unique(dat.lin.output)
  
#print out final data frame
  print(dat.lin.output)
}
