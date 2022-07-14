#
#
library(tidyverse)
library(readr)
###
#
#define input
dat.HQ.file <- "Intermediates/Analytical_Validation/AV_PPL_HQ_dat.csv"
Std.info.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
all.EE.file <- "Intermediates/Analytical_Validation/AV_PPL_EE_RF_Dat.csv"

##bring in and combine EE data
dat.HQ <- read_csv(dat.HQ.file) %>%
  rename("Compound.Name_old" = MF) %>%
  select(Compound.Name_old)

dat.EE <- read_csv(all.EE.file) %>%
  select(MF, Fraction, sample, Sample.Mean.EE, Sample.SD.EE) %>%
  unique() %>%
  rename("Compound.Name_old" = MF) %>%
  mutate(Compound.Name_old = str_replace_all(.$Compound.Name_old, "Isoleucine", "(Iso)leucine"))


###Join together, remove values in Puget Sound that are too high and replace with *ND
HQ.all.EE <- left_join(dat.HQ, dat.EE, by = "Compound.Name_old") %>%
  mutate(Sample.Mean.EE = case_when(Sample.Mean.EE >= 200 ~ "ND",
                              TRUE ~ as.character(Sample.Mean.EE))) %>% 
  mutate(Sample.SD.EE = case_when(Sample.SD.EE >= 200 ~ "ND",
                                    TRUE ~ as.character(Sample.SD.EE)))
##Bring in Standards info
Ing.name.dat <- read_csv(Std.info.file) %>%
  mutate(RT = RT..min.) %>%
  select(Compound.Name, Compound.Name_old, Compound.Name_figure, Column, z, Priority, m.z, RT, ionization_form) %>%
  mutate(Compound.Name_old = str_replace_all(.$Compound.Name_old, "Isoleucine", "(Iso)leucine"))


HQ.all.EE.2 <- left_join(HQ.all.EE, Ing.name.dat, by = "Compound.Name_old") %>%
  select(Compound.Name, Fraction, sample, Sample.Mean.EE, Sample.SD.EE) %>%
  unique()  %>%
  mutate(Compound.Name = str_replace_all(.$Compound.Name, "L-Isoleucine", "(Iso)leucine")) %>%
  mutate(Compound.Name = str_replace_all(.$Compound.Name, "3',5'-Cyclic AMP", "3',5'-Cyclic Adenosine Monophosphate")) %>%
  mutate(Compound.Name = str_replace_all(.$Compound.Name, "3',5'-Cyclic GMP", "3',5'-Cyclic Guanosine Monophosphate"))







###Rename and tidy up table for export
supp.table.5 <- HQ.all.EE.2 %>%
  mutate(Sample.Mean.EE = print(formatC(signif(as.numeric(Sample.Mean.EE), digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(Sample.SD.EE = print(formatC(signif(as.numeric(Sample.SD.EE), digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(Fraction = case_when(Fraction == "Pos" ~ "HILIC Pos",
                              Fraction == "Neg" ~ "HILIC Neg",
                              Fraction == "RP" ~ "RP")) %>%
  rename("Compound" = Compound.Name,
         'Sample' = sample,
         "Mean EE (%)" = Sample.Mean.EE,
         "Standard Deviation of EE (%)" = Sample.SD.EE
  ) %>%
  filter(!is.na(Compound))

#export
write_csv(supp.table.5, file = "Tables/Output/PPL_all_EE_supptable5.csv")
