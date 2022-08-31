#
#
library(tidyverse)
library(readr)
###
#
#define input
dat.HQ.file <- "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv"
Std.info.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
all.EE.file <- "Intermediates/Analytical_Validation/AV_CX_EE_RF_Dat.csv"
ANOVA.file <- "Intermediates/Analytical_Validation/AV_CX_ANOVA_Results.csv"

##bring in and combine EE data
dat.HQ <- read_csv(dat.HQ.file) %>%
  rename("Compound.Name_old" = MF) %>%
  select(Compound.Name_old, column)

dat.EE <- read_csv(all.EE.file) %>%
  select(MF, Fraction, sample, Sample.Mean.EE, Sample.SD.EE) %>%
  unique() %>%
  rename("Compound.Name_old" = MF) %>%
  mutate(Compound.Name_old = str_replace_all(.$Compound.Name_old, "Isoleucine", "(Iso)leucine"))

HQ.all.EE <- left_join(dat.HQ, dat.EE, by = "Compound.Name_old")


###ANOVA results
ANOVA.results <- read_csv(ANOVA.file) %>%
  select(MF, Pval.EE, Pval.rf) %>%
  unique() %>% 
  rename("Compound.Name_old" = MF)

HQ.all.EE.pvals <- left_join(HQ.all.EE, ANOVA.results)

##Bring in Standards info
Ing.name.dat <- read_csv(Std.info.file) %>%
  mutate(RT = RT..min.) %>%
  select(Compound.Name, Compound.Name_old, Compound.Name_figure, Column, z, Priority, m.z, RT, ionization_form)  %>%
  mutate(Compound.Name_old = str_replace_all(.$Compound.Name_old, "Isoleucine", "(Iso)leucine"))


HQ.all.EE.2 <- left_join(HQ.all.EE.pvals, Ing.name.dat, by = "Compound.Name_old") %>%
  select(Compound.Name, Fraction, sample, Sample.Mean.EE, Sample.SD.EE, Pval.EE, Pval.rf) %>%
  unique()  %>%
  mutate(Compound.Name = str_replace_all(.$Compound.Name, "L-Isoleucine", "(Iso)leucine")) %>%
  mutate(Compound.Name = str_replace_all(.$Compound.Name, "3',5'-Cyclic AMP", "3',5'-Cyclic Adenosine Monophosphate"))





###Rename and tidy up table for export
supp.table.4 <- HQ.all.EE.2 %>%
  mutate(sample = str_replace(.$sample, "A", "Aloha")) %>%
  mutate(Fraction = case_when(Fraction == "Pos" ~ "HILIC Pos",
                              Fraction == "Neg" ~ "HILIC Neg",
                              Fraction == "RP" ~ "RP")) %>%
  mutate(Sample.Mean.EE = print(formatC(signif(Sample.Mean.EE,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(Sample.SD.EE = print(formatC(signif(Sample.SD.EE,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(Pval.EE = print(formatC(signif(Pval.EE,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(Pval.rf = print(formatC(signif(Pval.rf,digits=3), digits=3,format="fg", flag="#"))) %>%
  rename("Compound" = Compound.Name,
         'Sample' = sample,
         "Mean EE (%)" = Sample.Mean.EE,
         "Standard Deviation of EE (%)" = Sample.SD.EE,
         "EE ANOVA Adjusted P-value" = Pval.EE,
         "RF ANOVA Adjusted P-value" = Pval.rf
         ) %>%
  arrange(Compound) %>%
  filter(!is.na(Compound))

#export
write_csv(supp.table.4, file = "Tables/Output/CX_all_EE_supptable4.csv")
