#THis script...
#
#
#
library(tidyverse)
library(readr)
###
#
#define input
#dat.EE.file <- "Intermediates/Analytical_Validation/AV_PPL_HQ_dat.csv"
#dat.HILIC.LOD.file <- "Intermediates/Environmental_Samples/ES_CX_Blk_LOD_Concentrations.csv"
#dat.RP.LOD.file <- "Intermediates/Environmental_Samples/ES_CXC_RP_Blk_LOD_Concentrations.csv"
#Std.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"

####Define inputs
HILIC.conc.file <- "Intermediates/Environmental_Samples/ES_PPL_HILIC_Concentrations.csv"
RP.conc.file <- "Intermediates/Environmental_Samples/ES_PPL_RP_Concentrations.csv"
HILIC.LOD.file <- "Intermediates/Environmental_Samples/ES_PPL_Blk_LOD_Concentrations.csv"
RP.LOD.file <- "Intermediates/Environmental_Samples/ES_PPL_RP_BLk_LOD_Concentrations.csv"
###
HILIC.conc <- read_csv(HILIC.conc.file) %>%
  mutate(sample = case_when(.$sample == "A" ~ "Aloha",
                            TRUE ~ sample))

RP.conc <- read_csv(RP.conc.file) 

##
dat <- rbind(HILIC.conc, RP.conc) %>%
  filter(!Compound == "Salicylic Acid") %>%
  group_by(Compound, sample) %>%
  mutate(mean.conc = mean(EE.adjust.conc),
         sd.mean.conc = sd(EE.adjust.conc),
         mean.Nmol.C = mean(Nmol.C),
         sd.Nmol.C = sd(Nmol.C),
         mean.Nmol.N = mean(Nmol.N),
         sd.Nmol.N = sd(Nmol.N)) %>%
  select(-EE.adjust.conc, -Nmol.C, -Nmol.N, -Rep) %>%
  unique()

####tidy up data for supplemental table
supp.dat <- dat %>%
  select(-C, -N) %>%
  rename("Sample" = sample,
         "Mean Concentration (nM)" = mean.conc,
         "Standard Deviation of Concentration (nM)" = sd.mean.conc,
         "Mean Concentration C (nM C)" = mean.Nmol.C,
         "Standard Deviation of Concentration C (nM C)" = sd.Nmol.C,
         "Mean Concentration N (nM N)" = mean.Nmol.N,
         "Standard Deviation of Concentration N (nM N)" = sd.Nmol.N)


###Bring in LOD values:
HILIC.LOD <- read_csv(HILIC.LOD.file) %>%
  select(Compound, EE.adjust.lod) %>%
  rename("Applied LOD threshold (nM)" = EE.adjust.lod)

RP.LOD <- read_csv(RP.LOD.file)  %>%
  select(Compound, EE.adjust.lod) %>%
  rename("Applied LOD threshold (nM)" = EE.adjust.lod)

All.LOD <- rbind(HILIC.LOD, RP.LOD)

###Add in Applied LOD values 
supp.dat.2 <- left_join(supp.dat, All.LOD)

#export:
write_csv(supp.dat.2, file = "Tables/Output/PPL_EnviroConc_supptable7.csv")
