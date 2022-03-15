#
#
#
#
#
library(tidyverse)
library(readr)


##
#Define inputs:
CX.file <- "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv"
PPL.file <- "Intermediates/Analytical_Validation/AV_PPL_HQ_dat.csv" 

###load in CX and PPL data and add method identifier 
CX.dat <- read_csv(CX.file) %>%
  mutate(method = "CX-SPE")
PPL.dat <- read_csv(PPL.file) %>%
  mutate(method = "PPL-SPE")

#combine data
comb.dat <- full_join(CX.dat, PPL.dat) %>%
  select(-R2, -z, -column)

###summarize total numbers of compounds extracted by each method in each fraction
sum.dat <- comb.dat %>%
  group_by(method, Fraction) %>%
  summarise(count = n())


###Examine level of overlap between CX-SPE and PPL-SPE
overlap.dat <- comb.dat %>%
  select(MF) %>%
  mutate(dup = duplicated(.$MF)) %>%
  group_by(dup) %>%
  summarize(count = n())

####Write comb.dat to a csv file
write_csv(comb.dat, path = "Intermediates/Analytical_Validation/AV_CX_PPL_HQ_combined.csv")





















