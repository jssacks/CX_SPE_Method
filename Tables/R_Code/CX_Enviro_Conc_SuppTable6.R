#
#
#
#
library(tidyverse)
library(readr)
#
#
####Define inputs
HILIC.conc.file <- "Intermediates/Environmental_Samples/ES_CX_HILIC_Concentrations.csv"
RP.conc.file <- "Intermediates/Environmental_Samples/ES_CXC_RP_Concentrations.csv"
#HILIC.LOD.file <- "Intermediates/Environmental_Samples/ES_Blk_LOD_Concentrations.csv"
#RP.LOD.file <- "Intermediates/Environmental_Samples/ES_CXC_RP_BLk_LOD_Concentrations.csv"
###
HILIC.conc <- read_csv(HILIC.conc.file) %>%
  mutate(sample = case_when(.$sample == "A" ~ "Aloha",
                            TRUE ~ sample))

RP.conc <- read_csv(RP.conc.file) %>%
  mutate(sample = str_replace_all(.$samp, "A", "Aloha")) %>%
  select(-samp)

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
         
#export:
write_csv(supp.dat, file = "Tables/Output/CX_EnviroConc_supptable6.csv")


