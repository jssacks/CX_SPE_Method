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
  select(Compound, Fraction, Overall.Mean.EE, Overall.RSD, R2, EE.adjust.lod) %>%
  mutate(Overall.Mean.EE = print(formatC(signif(Overall.Mean.EE,digits=3), digits=3,format="fg", flag = "#"))) %>%
  mutate(Overall.RSD = print(formatC(signif(Overall.RSD,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(EE.adjust.lod = print(formatC(signif(EE.adjust.lod,digits=2), digits=2,format="fg", flag="#"))) %>%
  mutate(R2 = print(formatC(signif(R2,digits=3), digits=3,format="fg", flag="#"))) %>%
  rename("Extraction Efficiency (%)" = Overall.Mean.EE,
         "RSD of EE (%)" = Overall.RSD,
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
  rename(Compound = Compound.Name) %>%
  mutate(Compound = str_replace_all(.$Compound, "L-Isoleucine", "(Iso)leucine")) %>%
  mutate(Compound = str_replace_all(.$Compound, "3',5'-Cyclic AMP", "3',5'-Cyclic Adenosine Monophosphate")) %>%
  arrange(Compound)
    
    
##write to a csv
write_csv(all.dat.std, file = "Tables/Output/Main_Text_1.csv")

