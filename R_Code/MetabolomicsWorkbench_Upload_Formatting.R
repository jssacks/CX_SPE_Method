#
#
#
#
#
#script for metabolomics workbench upload


library(tidyverse)
library(readr)
#
#Define inputs

##Standards Sheet
stds.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"

#HILIC
CX.HILIC.file <- "Intermediates/Environmental_Samples/ES_CatEx_HILIC_BMISed_BlkQCed_dat.csv"
PPL.HILIC.file <- "Intermediates/Environmental_Samples/ES_PPL_HILIC_BMISed_BlkQCed_dat.csv"

#RP





###Organize Standards info
stds.1 <- read_csv(stds.file) %>%
  select(Compound.Name_old, Emperical.Formula, ionization_form, Column, `RT..min.`, m.z, C0) %>%
  rename("Metabolite_name" = Compound.Name_old,
         "RT_min" = `RT..min.`,
         "KEGG_ID" = C0)



#HILIC Pos org
cx.hilic.dat <- read_csv(CX.HILIC.file)
ppl.hilic.dat <- read_csv(PPL.HILIC.file) 

hilic.pos.dat <- rbind(cx.hilic.dat, ppl.hilic.dat) %>%
  filter(!str_detect(.$MF, "_Neg")) %>%
  mutate(MF = str_replace_all(.$MF, "_Pos", "")) %>%
  mutate(SampName = str_replace_all(.$SampID, "210322_Smp_", "")) %>%
  mutate(SampName = str_replace_all(.$SampName, "210322_Poo_", "")) %>%
  mutate(SampName = str_replace_all(.$SampName, "Aloha", "ALOHA")) %>%
  select(MF, SampName, Adjusted_Area) %>%
  pivot_wider(id_cols = MF, names_from = SampName, values_from = Adjusted_Area)

write_tsv(hilic.pos.dat, file = "Intermediates/Environmental_Samples/CXSPE_MetabWorkbench_HILICPos.tsv")

###make hilic pos standards tsv
hilic.pos.stds <- stds.1 %>%
  filter(.$Metabolite_name %in% hilic.pos.dat$MF) %>%
  filter(!ionization_form == "[M-H]") %>%
  filter(Column == "HILIC")

write_tsv(hilic.pos.stds, file = "Intermediates/Environmental_Samples/CXSPE_MW_MetabMetadata_HILICPos.tsv")












