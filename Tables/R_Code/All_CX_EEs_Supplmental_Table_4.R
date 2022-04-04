#
#
library(tidyverse)
library(readr)
###
#
#define input
dat.HQ.file <- "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv"
Std.info.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
all.EE.file <- "Intermediates/Analytical_Validation/AV_CX_EE_RF_Dat.csv"

##bring in and combine EE data
dat.HQ <- read_csv(dat.HQ.file) %>%
  rename("Compound.Name_old" = MF) %>%
  select(Compound.Name_old, column)

dat.EE <- read_csv(all.EE.file) %>%
  select(MF, Fraction, sample, Sample.Mean.EE, Sample.SD.EE) %>%
  unique() %>%
  rename("Compound.Name_old" = MF)

HQ.all.EE <- left_join(dat.HQ, dat.EE, by = "Compound.Name_old")

##Bring in Standards info
Ing.name.dat <- read_csv(Std.info.file) %>%
  mutate(RT = RT..min.) %>%
  select(Compound.Name, Compound.Name_old, Compound.Name_figure, Column, z, Priority, m.z, RT, ionization_form)

HQ.all.EE.2 <- left_join(HQ.all.EE, Ing.name.dat, by = "Compound.Name_old") %>%
  select(Compound.Name, Fraction, sample, Sample.Mean.EE, Sample.SD.EE) %>%
  unique()



###Rename and tidy up table for export
supp.table.4 <- HQ.all.EE.2 %>%
  mutate(Fraction = case_when(Fraction == "Pos" ~ "HILIC Pos",
                              Fraction == "Neg" ~ "HILIC Neg",
                              Fraction == "RP" ~ "RP")) %>%
  rename("Compound" = Compound.Name,
         'Sample' = sample,
         "Mean EE (%)" = Sample.Mean.EE,
         "Standard Deviation of EE (%)" = Sample.SD.EE
         )

#export
write_csv(supp.table.4, file = "Tables/Output/CX_all_EE_supptable4.csv")
