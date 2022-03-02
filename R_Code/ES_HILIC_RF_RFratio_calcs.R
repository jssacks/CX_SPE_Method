###
#
#
##
#
#
#
#
#

###
library(readr)
library(tidyverse)

#Define Inputs:
HILIC.data.file <- "Intermediates/Environmental_Samples/ES_HILIC_targeted_combined_raw.csv"
Stds.info.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"

####Load in standards
HILIC_stds_dat <- read_csv(HILIC.data.file) %>%
  pivot_longer(cols = -Compound, names_to = "Rep", values_to = "Area") %>%
  filter(str_detect(.$Rep, "Std")) %>% 
  mutate(Compound = str_replace_all(.$Compound, "butyric_acid_Neg", "butyric acid_Neg")) %>%
  separate(Compound, into = c("Compound", "ion_mode"), sep = "_") %>%
  mutate(z = ifelse(ion_mode == "Neg", -1, 1)) %>%
  select(Compound, Rep, z, Area) %>%
  mutate(Mix = str_extract(Rep, "Mix\\d"))

###Get Mix and Concentration info:
Stds.info <- read_csv(Stds.info.file) %>%
  filter(Priority == TRUE) %>%
  select(Column, Compound.Name_old, z, HILICMix, Conc..uM) %>%
  rename(Compound = Compound.Name_old,
         Conc = Conc..uM) %>%
  mutate(Compound = str_replace_all(.$Compound, "butyric_acid_Neg", "butyric acid_Neg"))

HILIC.stds.info <- Stds.info %>%
  filter(Column == "HILIC") %>%
  select(Compound, z, HILICMix, Conc)

###Join stuff together + remove Matrix Samples
HILIC.stds.dat.info <- full_join(HILIC_stds_dat, HILIC.stds.info, by = c("Compound", "z")) 

##Calculate RFs
HILIC.RFs <- HILIC.stds.dat.info %>%
  filter(Mix == HILICMix) %>%
  select(-Mix, -HILICMix) %>%
  filter(!str_detect(.$Rep, "Matrix")) %>%
  mutate(RF = as.numeric(Area)/Conc, NA) %>%
  group_by(Compound, z) %>%
  summarise(RFmax = max(RF),
            RFmin = min(RF),
            RF = mean(RF, na.rm = TRUE))

HILIC.RFratio.dat <- HILIC.stds.dat.info %>%
  filter(HILICMix == Mix | is.na(Mix)) %>% 
  mutate(RunNumber = str_extract(Rep, "_\\d$")) %>%
  mutate(RunType = ifelse(str_detect(Rep, "StdsMix\\dInH2O")|
                            str_detect(Rep, "StdsInH2O"), "Std_in_h2O", 
                          ifelse(str_detect(Rep, "StdsMix\\dInMatrix") |
                                   str_detect(Rep, "StdsInMatrix"), "Std_in_matrix",
                                 "Matrix_in_h2O"))) %>%
  filter(HILICMix == Mix | is.na(Mix)) %>% 
  select(-Mix, -HILICMix, -Rep, -Conc) %>%
  spread(key = RunType, value = Area ) 

HILIC.RF.ratio <- HILIC.RFratio.dat %>%
  ungroup() %>%
  group_by(Compound, RunNumber, z) %>%
  summarize("Std_in_matrix" = Std_in_matrix,
            "Matrix_in_h2o" = Matrix_in_h2O,
            "Std_in_h2o" = Std_in_h2O) %>%
  mutate(RFratio = (Std_in_matrix - Matrix_in_h2o)/ Std_in_h2o) %>%
  group_by(Compound, z) %>%
  summarise(RFratio = mean(RFratio))

###Join it all together
HILIC.dat.RF.RFratio <- left_join(HILIC.RFs, HILIC.RF.ratio)
write_csv(HILIC.dat.RF.RFratio, file = "Intermediates/Environmental_Samples/HILIC_RFs_RFratios.csv")
