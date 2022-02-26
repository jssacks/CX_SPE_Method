
#
# This Script takes in BMIS normalized data from RP and HILIC and combines the data into one
# data frame while removing duplicate compounds in the data by only selecting for the optimal column + polarity
# for each compound. The output files are combined HILIC and RP data files for both CX and PPL experiments
#
#
#
#
library(tidyverse)
library(readr)
###


#Define inputs:
CX.HILIC.dat <- "Intermediates/Analytical_Validation/AV_CX_HILIC_BMISed_dat.csv"
CX.RP.dat <- "Intermediates/Analytical_Validation/AV_CX_RP_BMISed_dat.csv"
PPL.HILIC.dat <- "Intermediates/Analytical_Validation/AV_PPL_HILIC_BMISed_dat.csv"
PPL.RP.dat <- "Intermediates/Analytical_Validation/AV_PPL_RP_BMISed_dat.csv"
Std.info.dat <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
Std.Spike.dat <- "Meta_Data/Ingalls_Standards/Ingalls_Standards_AV_SpikeConc.csv"
####

###Clean up and organize and combine CX HILIC and RP Data while adding a column designation 
  # to distinguish the two data sets
dat.HILIC <- read_csv(CX.HILIC.dat) %>%
  mutate(MF = str_replace_all(MF, pattern="yric_acid", 
                              replacement = "yric acid")) %>%
  separate(MF, into = c("MF", "Fraction"), sep = "_") %>%
  mutate("column" = "HILIC")

dat.RP <- read_csv(CX.RP.dat) %>%
  mutate("column" = "RP") %>%
  mutate("Fraction" = "RP") 

dat.all <- rbind(dat.HILIC, dat.RP) %>%
  mutate(z = if_else(Fraction == "Neg", -1, 1))


###Load in Ingalls lab standard and AV spike concentration data, clean up and remove extraneous rows and columns
spike.conc.dat <- read_csv(Std.Spike.dat) %>%
  filter(!str_detect(.$Compound_Name, "skyline")) %>%
  rename("Compound.Name" = Compound_Name)
Ing.name.dat <- read_csv(Std.info.dat) %>%
  select(Compound.Name, Compound.Name_old, Column, z, Priority)

stds.dat <- left_join(Ing.name.dat, spike.conc.dat) %>%
  filter(Priority == TRUE) %>%
  drop_na() %>%
  mutate("MF" = Compound.Name_old) %>%
  mutate("column" = Column) %>%
  select(MF, column, z, Priority, Spike_Concentration_uM)
  


###Deduplicate data by only selecting the data for the each compound from its optimal column and polarity
dat.dedup <- left_join(dat.all, stds.dat, by = c("MF", "column", "z")) %>%
filter(Priority == TRUE) 
write_csv(dat.dedup, path = "Intermediates/Analytical_Validation/AV_CX_Full_Dat.csv")








# PPL HILIC and RP Combine ------------------------------------------------

###Clean up and organize and combine PPL HILIC and RP Data while adding a column designation 
# to distinguish the two data sets
PPL.dat.HILIC <- read_csv(PPL.HILIC.dat) %>%
  mutate(MF = str_replace_all(MF, pattern="yric_acid", 
                              replacement = "yric acid")) %>%
  separate(MF, into = c("MF", "Fraction"), sep = "_") %>%
  mutate("column" = "HILIC")

PPL.dat.RP <- read_csv(PPL.RP.dat) %>%
  mutate_if(.,is.character,str_replace_all, pattern="PPL_", replacement = "") %>%
  mutate("column" = "RP") %>%
  mutate("Fraction" = "RP") 

PPL.dat.all <- rbind(PPL.dat.HILIC, PPL.dat.RP) %>%
  mutate(z = if_else(Fraction == "Neg", -1, 1))


###Load in Ingalls lab standard and AV spike concentration data, clean up and remove extraneous rows and columns
spike.conc.dat <- read_csv(Std.Spike.dat) %>%
  filter(!str_detect(.$Compound_Name, "skyline")) %>%
  rename("Compound.Name" = Compound_Name)
Ing.name.dat <- read_csv(Std.info.dat) %>%
  select(Compound.Name, Compound.Name_old, Column, z, Priority)

stds.dat <- left_join(Ing.name.dat, spike.conc.dat) %>%
  filter(Priority == TRUE) %>%
  drop_na() %>%
  mutate("MF" = Compound.Name_old) %>%
  mutate("column" = Column) %>%
  select(MF, column, z, Priority, Spike_Concentration_uM)

###Deduplicate data by only selecting the data for the each compound from its optimal column and polarity
PPL.dat.dedup <- left_join(PPL.dat.all, stds.dat, by = c("MF", "column", "z")) %>%
  filter(Priority == TRUE) %>%
  mutate_if(.,is.character,str_replace_all, pattern="SBf_", replacement = "SBe_")
write_csv(PPL.dat.dedup, path = "Intermediates/Analytical_Validation/AV_PPL_Full_Dat.csv")
  
  
  
  
  
  
  
























