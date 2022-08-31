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
Std.info.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
###
HILIC.conc <- read_csv(HILIC.conc.file) %>%
  mutate(sample = case_when(.$sample == "A" ~ "Aloha",
                            TRUE ~ sample))

RP.conc <- read_csv(RP.conc.file) 

##cleaning up data and finding means. also removing salicylic acid and L-pyroglutamic acid because the data is unreliable
conc.dat <- rbind(HILIC.conc, RP.conc) %>%
  filter(!Compound == "Salicylic Acid") %>%
  filter(!Compound == "L-Pyroglutamic acid") %>%
  group_by(Compound, sample) %>%
  mutate(mean.conc = mean(EE.adjust.conc),
         sd.mean.conc = sd(EE.adjust.conc),
         mean.Nmol.C = mean(Nmol.C),
         sd.Nmol.C = sd(Nmol.C),
         mean.Nmol.N = mean(Nmol.N),
         sd.Nmol.N = sd(Nmol.N)) %>%
  select(-EE.adjust.conc, -Nmol.C, -Nmol.N, -Rep) %>%
  unique() %>%
  rename("Compound.Name_old" = Compound) %>%
  filter(!sample == "Mort")


###Fix Names:
##Bring in Standards info
Ing.name.dat <- read_csv(Std.info.file) %>%
  mutate(RT = RT..min.) %>%
  select(Compound.Name, Compound.Name_old)  %>%
  mutate(Compound.Name_old = str_replace_all(.$Compound.Name_old, "Isoleucine", "(Iso)leucine")) %>%
  unique()


conc.dat.2 <- left_join(conc.dat, Ing.name.dat, by = "Compound.Name_old") %>%
  ungroup() %>%
  mutate(Compound.Name = str_replace_all(.$Compound.Name, "L-Isoleucine", "(Iso)leucine")) %>%
  mutate(Compound.Name = str_replace_all(.$Compound.Name, "3',5'-Cyclic AMP", "3',5'-Cyclic Adenosine Monophosphate")) %>%
  rename("Compound" = Compound.Name) %>%
  select(-Compound.Name_old) %>%
  select(Compound, everything())






####tidy up data for supplemental table
supp.dat <- conc.dat.2%>%
  select(-C, -N) %>%
  mutate(mean.conc = print(formatC(signif(mean.conc,digits=2), digits=2,format="fg", flag="#"))) %>%
  mutate(sd.mean.conc = print(formatC(signif(sd.mean.conc,digits=2), digits=2,format="fg", flag="#"))) %>%
  mutate(mean.Nmol.C = print(formatC(signif(mean.Nmol.C,digits=2), digits=2,format="fg", flag="#"))) %>%
  mutate(sd.Nmol.C = print(formatC(signif(sd.Nmol.C,digits=2), digits=2,format="fg", flag="#"))) %>%
  mutate(mean.Nmol.N= print(formatC(signif(mean.Nmol.N,digits=2), digits=2,format="fg", flag="#"))) %>%
  mutate(sd.Nmol.N= print(formatC(signif(sd.Nmol.N,digits=2), digits=2,format="fg", flag="#"))) %>%
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
  mutate(EE.adjust.lod = print(formatC(signif(EE.adjust.lod,digits=2), digits=2,format="fg", flag="#"))) %>%
  rename("Applied LOD threshold (nM)" = EE.adjust.lod)

RP.LOD <- read_csv(RP.LOD.file)  %>%
  select(Compound, EE.adjust.lod) %>%
  mutate(EE.adjust.lod = print(formatC(signif(EE.adjust.lod,digits=2), digits=2,format="fg", flag="#"))) %>%
  rename("Applied LOD threshold (nM)" = EE.adjust.lod)

All.LOD <- rbind(HILIC.LOD, RP.LOD) %>%
  rename("Compound.Name_old" = Compound) 

All.LOD.name <- left_join(All.LOD, Ing.name.dat, by = "Compound.Name_old") %>%
                            ungroup() %>%
                            mutate(Compound.Name = str_replace_all(.$Compound.Name, "L-Isoleucine", "(Iso)leucine")) %>%
                            mutate(Compound.Name = str_replace_all(.$Compound.Name, "3',5'-Cyclic AMP", "3',5'-Cyclic Adenosine Monophosphate")) %>%
                            rename("Compound" = Compound.Name) %>%
                            select(-Compound.Name_old) %>%
                            select(Compound, everything())

###Add in Applied LOD values 
supp.dat.2 <- left_join(supp.dat, All.LOD.name) %>%
  arrange(Compound)

#export:
write_csv(supp.dat.2, file = "Tables/Output/PPL_EnviroConc_supptable7.csv")

