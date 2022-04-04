



#
#
#
library(tidyverse)
library(readr)


Std.info.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
Std.Spike.file <- "Meta_Data/Ingalls_Standards/Ingalls_Standards_AV_SpikeConc.csv"
Supplier.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_Suppliers.csv"


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

####Change names and select final columns and column order to include in supplemental Table 1
supp.table.1 <- stds.dat.2 %>%
  select(Compound.Name, m.z, ionization_form, Column, RT, Spike_Concentration_uM, Supplier, Part_Number) %>%
  rename("Compound" = Compound.Name,
         "m/z" = m.z,
         "Ionization form" = ionization_form,
         "Retention Time" = RT,
         "Spike Concentration (uM)" = Spike_Concentration_uM,
         "Part Number" = Part_Number)
###Export
write_csv(supp.table.1, file = "Tables/Output/Std_Details_supptable1.csv")




