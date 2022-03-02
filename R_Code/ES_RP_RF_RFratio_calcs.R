###
library(readr)
library(tidyverse)

#Define Inputs:
RP.data.file <- "Intermediates/Environmental_Samples/ES_CX_RP_targeted_combined_raw.csv"
Stds.info.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"

####Load in standards
RP_stds_dat <- read_csv(RP.data.file) %>%
  pivot_longer(cols = -Compound, names_to = "Rep", values_to = "Area") %>%
  filter(str_detect(.$Rep, "Std")) %>% 
  mutate(Compound = str_replace_all(.$Compound, "butyric_acid_Neg", "butyric acid_Neg")) %>%
  separate(Compound, into = c("Compound", "ion_mode"), sep = "_") %>%
  mutate(z = ifelse(ion_mode == "Neg", -1, 1)) %>%
  select(Compound, Rep, Area) 

###Make areas that are NAs equal to 0
RP_stds_dat$Area[is.na(RP_stds_dat$Area)] <- 0

###Get Mix and Concentration info:
Stds.info <- read_csv(Stds.info.file) %>%
  filter(Priority == TRUE) %>%
  select(Column, Compound.Name_old, z, HILICMix, Conc..uM) %>%
  rename(Compound = Compound.Name_old,
         Conc = Conc..uM) %>%
  mutate(Compound = str_replace_all(.$Compound, "butyric_acid_Neg", "butyric acid_Neg"))

RP.stds.info <- Stds.info %>%
  filter(Column == "RP") %>%
  select(Compound, z, Conc)

###Join stuff together + remove Matrix Samples
RP.stds.dat.info <- full_join(RP_stds_dat, RP.stds.info, by = c("Compound")) 

##Calculate RFs
RP.RFs <- RP.stds.dat.info %>%
  filter(!str_detect(.$Rep, "Matrix")) %>%
  mutate(RF = as.numeric(Area)/Conc, NA) %>%
  group_by(Compound, z) %>%
  summarise(RFmax = max(RF),
            RFmin = min(RF),
            RF = mean(RF, na.rm = TRUE))

RP.RFratio.dat <- RP.stds.dat.info %>%
  #  filter(HILICMix == Mix | is.na(Mix)) %>% 
  mutate(RunNumber = str_extract(Rep, "_\\d$")) %>%
  mutate(RunType = ifelse(str_detect(Rep, "StdsMix\\dInH2O")|
                            str_detect(Rep, "StdsInH2O"), "Std_in_h2O", 
                          ifelse(str_detect(Rep, "StdsMix\\dInMatrix") |
                                   str_detect(Rep, "StdsInMatrix"), "Std_in_matrix",
                                 "Matrix_in_h2O"))) %>%
  # filter(HILICMix == Mix | is.na(Mix)) %>% 
  select(-Rep, -Conc) %>%
  spread(key = RunType, value = Area) 

RP.RF.ratio <- RP.RFratio.dat %>%
  ungroup() %>%
  group_by(Compound, RunNumber, z) %>%
  summarize("Std_in_matrix" = Std_in_matrix,
            "Matrix_in_h2o" = Matrix_in_h2O,
            "Std_in_h2o" = Std_in_h2O) %>%
  mutate(RFratio = (Std_in_matrix - Matrix_in_h2o)/ Std_in_h2o) %>%
  group_by(Compound, z) %>%
  summarise(RFratio = mean(RFratio))

###Join it all together
RP.dat.RF.RFratio <- left_join(RP.RFs, RP.RF.ratio)
write_csv(RP.dat.RF.RFratio, path = "Intermediates/Environmental_Samples/ES_CX_RP_RFs_RFratios.csv")
