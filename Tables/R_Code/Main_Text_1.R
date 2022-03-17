#
#
#
#THis script...
#
#
#
library(tidyverse)
library(readr)
###
#
#define input
dat.EE.file <- "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv"
dat.HILIC.LOD.file <- "Intermediates/Environmental_Samples/ES_CX_Blk_LOD_Concentrations.csv"
dat.RP.LOD.file <- "Intermediates/Environmental_Samples/ES_CXC_RP_Blk_LOD_Concentrations.csv"
Std.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"

##read in EE data
dat.EE <- read_csv(dat.EE.file) %>%
  rename("Compound" = MF)

##read in LOD data
dat.HILIC.LOD <- read_csv(dat.HILIC.LOD.file) %>%
  select(-EE.adjust.loq)
dat.RP.LOD <- read_csv(dat.RP.LOD.file)

#combine LOD data
all.LOD.dat <- rbind(dat.HILIC.LOD, dat.RP.LOD)

###combine EE and LOD data, select only desired columns, tidy up, rename things to be appropriate 
all.dat <- left_join(dat.EE, all.LOD.dat) %>%
  select(Compound, Fraction, Sample.Mean.EE, Sample.RSD, R2, EE.adjust.lod) %>%
  rename("Extraction Efficiency (%)" = Sample.Mean.EE,
         "RSD of EE (%)" = Sample.RSD,
         "LOD (nM)" = EE.adjust.lod) %>%
  mutate(Fraction = case_when(
    .$Fraction == "Pos" ~ "HILIC Pos",
    .$Fraction == "Neg" ~ "HILIC Neg",
    TRUE ~ Fraction
  )) %>%
  arrange(Compound)

###Bring in Standards Information, get just new compound name information
std.info <- read_csv(Std.file)
std.info.2 <- std.info %>%
  filter(Priority = TRUE) %>%
  rename(Compound = Compound.Name_old) %>%
  mutate(Compound = str_replace_all(.$Compound, "Isoleucine", "(Iso)leucine")) %>%
  filter(!Compound == "Leucine") %>%
  select(Compound, Compound.Name)

###Combine standards information with EE+LOD Table
all.dat.std <- left_join(all.dat, std.info.2) %>%
  unique() %>%
  select(Compound.Name, Fraction, `Extraction Efficiency (%)`, `RSD of EE (%)`, R2, `LOD (nM)`) %>%
  rename(Compound = Compound.Name)

##write to a csv
write_csv(all.dat.std, file = "Tables/Output/Main_Text_1.csv")
