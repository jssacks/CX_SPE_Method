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
Std.info.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
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
conc.dat <- rbind(HILIC.conc, RP.conc) %>%
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
  rename("Compound.Name_old" = Compound)

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







####tidy up data for supplemental table including sig figs
supp.dat <- conc.dat.2 %>%
  mutate(mean.conc = print(formatC(signif(mean.conc,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(sd.mean.conc = print(formatC(signif(sd.mean.conc,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(mean.Nmol.C = print(formatC(signif(mean.Nmol.C,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(sd.Nmol.C = print(formatC(signif(sd.Nmol.C,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(mean.Nmol.N= print(formatC(signif(mean.Nmol.N,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(sd.Nmol.N= print(formatC(signif(sd.Nmol.N,digits=3), digits=3,format="fg", flag="#"))) %>%
  select(-C, -N) %>%
  mutate(sample = str_replace(.$sample, "Alohaloha", "Aloha")) %>%
  rename("Sample" = sample,
         "Mean Concentration (nM)" = mean.conc,
         "Standard Deviation of Concentration (nM)" = sd.mean.conc,
         "Mean Concentration C (nM C)" = mean.Nmol.C,
         "Standard Deviation of Concentration C (nM C)" = sd.Nmol.C,
         "Mean Concentration N (nM N)" = mean.Nmol.N,
         "Standard Deviation of Concentration N (nM N)" = sd.Nmol.N) %>%
  arrange(Compound)
         
#export:
write_csv(supp.dat, file = "Tables/Output/CX_EnviroConc_supptable6.csv")


