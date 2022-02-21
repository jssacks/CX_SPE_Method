#####
##This script imports the raw skyline output files for the HILIC CX-SPE analytical validation (AV)
# section of this project, cleans, renames columns, and organizes the data in preparation
# for normalization. This script also cleans up the sample list for this data set and makes sure it 
# is consistent with the skyline output

#load packages
library(tidyverse)
library(readr)
library(stringr)

#load functions
source("R_Code/CXSPE_Functions.R")

#Define inputs
#HILIC positive skyline AV files for CX-SPE
g1.pos <- "Raw_Data/Skyline_Output/Analytical_Validation/CXSPE_AV_RP_G1.csv"
g2.pos <- "Raw_Data/Skyline_Output/Analytical_Validation/CXSPE_AV_RP_G2.csv"
g3.pos <- "Raw_Data/Skyline_Output/Analytical_Validation/CXSPE_AV_RP_G3.csv"
g4.pos <- "Raw_Data/Skyline_Output/Analytical_Validation/CXSPE_AV_RP_G4.csv"

#sample list
smp.list <- "Meta_Data/Sample_Lists/CXSPE_AV_RP_SampList.csv"


#################### ___________________________________________
#####load and parse files
g1.p <- csv_parse(g1.pos)
g2.p <- csv_parse(g2.pos)
g3.p <- csv_parse(g3.pos)
g4.p <- csv_parse(g4.pos)

#Need to add "Smp_" and "Poo_" to g4.p to fix naming issues
g4.p.LinA <- g4.p %>%
  filter(str_detect(g4.p$Rep, "LinA_0")) %>%
  mutate(Smp = "Smp") %>%
  unite("Rep", c("Smp", "Rep"))

g4.p.Poo <- g4.p %>%
  filter(str_detect(g4.p$Rep, "Poo")) %>%
  mutate(Poo = "Poo") %>%
  unite("Rep", c("Poo", "Rep"))

#Combine groups together
full.rp <- bind_rows(g1.p, g2.p, g3.p, g4.p.LinA, g4.p.Poo)

##remove blanks and standard mixes, add in date
full.rp.ns <- full.rp %>%
  filter(!str_detect(.$Rep, "Std_"))  %>%
  filter(!str_detect(.$Rep, "blk")) %>%
  filter(!str_detect(.$Rep, pattern = "Blk")) %>%
  filter(!str_detect(.$Rep, pattern = "PPL")) %>%
  mutate(date = "200825") %>%
  unite("Rep", c("date", "Rep")) %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point1_S_", replacement = "LinA_0point1_") %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point5_S_", replacement = "LinA_0point5_")

###remove internal standards, 
AV.RP.full <- full.rp.ns %>%
  filter(!str_detect(.$Compound, ","))

#pivot dataframe wider and write to a csv for normalization later 
AV.RP.full.export <- pivot_wider(AV.RP.full, names_from = Rep, values_from = Area)
write_csv(AV.RP.full.export, path = "Intermediates/Analytical_Validation/AV_RP_combined_raw.csv")

####Create IS dataframe
full.IS <- full.rp.ns %>%
  filter(str_detect(.$Compound, ",")) %>%
  filter(!Compound == "2,4 decadienal") %>%
  filter(!Compound == "2,4 octadienal")
write_csv(full.IS, path = "Intermediates/Analytical_Validation/AV_RP_IS.csv")


#####Make sample list compatible with data frame by removing blanks, PPL samples, adding dates and fixing names
rp.smp.list.AV <- read_csv(smp.list, skip = 1) %>%
  rename(Rep = `File Name`) %>%
  filter(!str_detect(.$Rep, "Blk")) %>%
  filter(!str_detect(.$Rep, "blk")) %>%
  filter(!str_detect(.$Rep, "PPL")) %>%
  mutate_if(.,is.character,str_replace_all, pattern="200826_", replacement = "") %>%
  mutate_if(.,is.character,str_replace_all, pattern="200812_", replacement = "") %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point1_S_", replacement = "LinA_0point1_") %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point5_S_", replacement = "LinA_0point5_") %>%
  filter(!str_detect(.$Rep, "NoSpike"))%>%
  mutate(date = "200825") %>%
  unite("Rep", c("date", "Rep")) %>%
  rename("Injec_vol" = InjVol)

###write sample list to  a data frame 
write_csv(rp.smp.list.AV, path = "Intermediates/Analytical_Validation/AV_CX_RP_Samp_List.csv")










