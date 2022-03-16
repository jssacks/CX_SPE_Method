#
#
#
#
#
library(readr)
library(tidyverse)
##
#Define inputs
HILIC.conc.file <- "Intermediates/Environmental_Samples/ES_CX_HILIC_Concentrations.csv"
HILIC.LOD.file <- "Intermediates/Environmental_Samples/ES_CX_Blk_LOD_Concentrations.csv"
RP.conc.file <- "Intermediates/Environmental_Samples/ES_CXC_RP_Concentrations.csv"
RP.LOD.file <- "Intermediates/Environmental_Samples/ES_CXC_RP_Blk_LOD_Concentrations.csv"
#
#
# Load in LOD data and combine 
HILIC.lod <- read_csv(HILIC.LOD.file) %>%
  select(-EE.adjust.loq)
RP.lod <- read_csv(RP.LOD.file)
all.lod <- rbind(HILIC.lod, RP.lod)


####Load in Environmental Concentrations 
HILIC.conc <- read_csv(HILIC.conc.file) %>%
  mutate(sample = case_when(.$sample == "A" ~ "Aloha",
                            TRUE ~ sample))

RP.conc <- read_csv(RP.conc.file) %>%
  mutate(sample = str_replace_all(.$samp, "A", "Aloha")) %>%
  select(-samp)

all.conc <- rbind(HILIC.conc, RP.conc) 


###summarize concentration data, #Salicyclic acid removed because of high and variable blanks
Conc.means <- all.conc%>%
  filter(!Compound == "Salicylic Acid") %>%
  group_by(Compound, sample) %>%
  mutate(mean.conc = mean(EE.adjust.conc),
         sd.mean.conc = sd(EE.adjust.conc),
         mean.Nmol.C = mean(Nmol.C),
         sd.Nmol.C = sd(Nmol.C),
         mean.Nmol.N = mean(Nmol.N),
         sd.Nmol.N = sd(Nmol.N)) %>%
  select(-EE.adjust.conc, -Nmol.C, -Nmol.N, -Rep) %>%
  unique() %>%
  ungroup()

Conc.sum <- Conc.means %>%
  group_by(sample) %>%
  summarise(count = n(),
            Mol.C = sum(mean.Nmol.C),
            Mol.N = sum(mean.Nmol.N, na.rm = TRUE))

