######
##
##This Script...
##
##

#packages
library(tidyverse)
library(readr)


##Define Inputs
CX.file <- "Intermediates/Analytical_Validation/AV_CX_Full_Dat.csv"
PPL.file <- "Intermediates/Analytical_Validation/AV_PPL_Full_Dat.csv"


####



#________CX-SPE EE and RF Calculations ---------------------------------------------------
CX.dat <- read_csv(CX.file) %>%
  select(MF, Fraction, z, column, SampID, Adjusted_Area, Spike_Concentration_uM) %>%
  separate(SampID, 
           c("runDate",
             "type","sample", "treatment", "replicate"),"_", remove = TRUE) %>%
  filter(type == "Smp") %>%
  filter(!sample == "LinA") %>%
  filter(!sample == "C") %>%
  unique()

###Calculate extraction efficiency values (EE) and response factor values in peak area/nM
CX.calc <- CX.dat %>%
  pivot_wider(id_cols = c("MF", "Fraction", "type", "sample", "replicate", "Spike_Concentration_uM"), names_from = treatment, values_from = Adjusted_Area) %>%
  mutate(EE = (SBe-NS)/(SAf-NS)*100,
         RF = (SAf-NS)/(Spike_Concentration_uM*1000))

##summarize EE and RF values
CX.sum <- CX.calc %>%
  group_by(MF) %>%
  mutate(Overall.Mean.EE = mean(EE),
         Overall.SD.EE = sd(EE),
         Overall.Mean.RF = mean(RF),
         Overall.SD.RF = sd(RF)) %>%
  group_by(MF, sample) %>%
  mutate(Sample.Mean.EE = mean(EE),
         Sample.SD.EE = sd(EE),
         Sample.Mean.RF = mean(RF),
         Sample.SD.RF = sd(RF))

####
write_csv(CX.sum, path = "Intermediates/Analytical_Validation/AV_CX_EE_RF_Dat.csv")



#________PPL-SPE EE and RF Calculations ---------------------------------------------------
PPL.dat <- read_csv(PPL.file) %>%
  select(MF, Fraction, z, column, SampID, Adjusted_Area, Spike_Concentration_uM) %>%
  separate(SampID, 
           c("runDate",
             "type","sample", "treatment", "replicate"),"_", remove = TRUE) %>%
  filter(type == "Smp") %>%
  filter(!sample == "LinA") %>%
  filter(!sample == "C") %>%
  unique()

###Calculate extraction efficiency values (EE) and response factor values in peak area/nM
PPL.calc <- PPL.dat %>%
  pivot_wider(id_cols = c("MF", "Fraction", "type", "sample", "replicate", "Spike_Concentration_uM"), names_from = treatment, values_from = Adjusted_Area) %>%
  mutate(EE = (SBe-NS)/(SAf-NS)*100,
         RF = (SAf-NS)/(Spike_Concentration_uM*1000))

##summarize EE and RF values
PPL.sum <- PPL.calc %>%
  group_by(MF) %>%
  mutate(Overall.Mean.EE = mean(EE),
         Overall.SD.EE = sd(EE),
         Overall.Mean.RF = mean(RF),
         Overall.SD.RF = sd(RF)) %>%
  group_by(MF, sample) %>%
  mutate(Sample.Mean.EE = mean(EE),
         Sample.SD.EE = sd(EE),
         Sample.Mean.RF = mean(RF),
         Sample.SD.RF = sd(RF))

####
write_csv(PPL.sum, path = "Intermediates/Analytical_Validation/AV_PPL_EE_RF_Dat.csv")


