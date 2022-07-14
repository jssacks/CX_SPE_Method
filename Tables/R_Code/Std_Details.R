



#
#
#
library(tidyverse)
library(readr)


Std.info.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
Std.Spike.file <- "Meta_Data/Ingalls_Standards/Ingalls_Standards_AV_SpikeConc.csv"
Supplier.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_Suppliers.csv"
CX.comps.file <- "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv"
PPL.comps.file <- "Intermediates/Analytical_Validation/AV_PPL_HQ_dat.csv"

###Load in Ingalls lab standard and AV spike concentration data, clean up and remove extraneous rows and columns
spike.conc.dat <- read_csv(Std.Spike.file) %>%
  filter(!str_detect(.$Compound_Name, "skyline")) %>%
  rename("Compound.Name" = Compound_Name)
Ing.name.dat <- read_csv(Std.info.file) %>%
  mutate(RT = RT..min.) %>%
  select(Compound.Name, Compound.Name_old, Compound.Name_figure, Column, z, Priority, m.z, RT, ionization_form)

###
stds.dat <- left_join(Ing.name.dat, spike.conc.dat) %>%
  filter(Priority == TRUE) %>%
  drop_na() %>%
  mutate("MF" = Compound.Name_old) %>%
  mutate("column" = Column) 


####Add in supplier information
supply.dat <- read_csv(Supplier.file) %>%
  mutate(MF = Compound_Name_Original) %>%
  select(MF, Supplier, Part_Number) %>%
  unique()

stds.dat.2 <- left_join(stds.dat, supply.dat)


#######Summarize details to include in manuscript

###summarize info by column
stds.summary <- stds.dat.2 %>%
  group_by(Column) %>%
  summarise(count = n(),
            med.spike = median(Spike_Concentration_uM)) 

###Add in info related to CX-SPE and PPL-SPE performance
CX.comps.dat <- read_csv(CX.comps.file) %>%
  select(MF) %>%
  mutate(CX = TRUE)
PPL.comps.dat <- read_csv(PPL.comps.file) %>%
  select(MF) %>%
  mutate(PPL = TRUE)

status.dat <- full_join(CX.comps.dat, PPL.comps.dat) %>%
  mutate(Status = case_when(
    CX == TRUE & PPL == TRUE ~ "CX-SPE and PPL-SPE",
    CX == TRUE & is.na(PPL) ~ "CX-SPE only",
    is.na(CX) & PPL == TRUE ~ "PPL-SPE only"))
status.dat.2 <- status.dat %>%
  select(MF, Status)

stds.dat.3 <- left_join(stds.dat.2, status.dat.2) %>%
  mutate(Status = case_when(
    is.na(Status) ~ "Neither Method",
    TRUE ~ Status
  ))

####Change names, change uM to nM, fix sig figs, and select final columns and column order to include in supplemental Table 1
supp.table.1 <- stds.dat.3 %>%
  mutate(Spike_Concentration_nM =1000*Spike_Concentration_uM) %>% 
  mutate(spike_conc_nM = print(formatC(signif(Spike_Concentration_nM,digits=3), digits=3,format="fg", flag="#"))) %>%
  mutate(RT = print(formatC(signif(RT,digits=3), digits=3,format="fg", flag="#"))) %>%
  select(Compound.Name, m.z, ionization_form, Column, RT, spike_conc_nM, Supplier, Part_Number, Status) %>%
  rename("Compound" = Compound.Name,
         "m/z" = m.z,
         "Ionization form" = ionization_form,
         "Retention Time" = RT,
         "Spike Concentration (nM)" = spike_conc_nM,
         "Part Number" = Part_Number) %>%
  unique() %>%
  arrange(Compound)
  
  
###Export
write_csv(supp.table.1, file = "Tables/Output/Std_Details_supptable1.csv")




