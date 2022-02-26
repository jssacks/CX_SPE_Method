#####
##This script imports the raw skyline output files for the HILIC PPL analytical validation (AV)
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
#HILIC positive skyline AV files for PPL
g1.pos <- "Raw_Data/Skyline_Output/Analytical_Validation/CXSPE_AV_Pos_G1.csv"
g6.pos <- "Raw_Data/Skyline_Output/Analytical_Validation/CXSPE_AV_Pos_G6.csv"


#HILIC negative skyline AV files for PPL
g1.neg <- "Raw_Data/Skyline_Output/Analytical_Validation/CXSPE_AV_Neg_G1.csv"
g6.neg <- "Raw_Data/Skyline_Output/Analytical_Validation/CXSPE_AV_Neg_G6.csv"

#
smp.list <- "Meta_Data/Sample_Lists/CXSPE_AV_HILIC_SampList.csv"


#################### ___________________________________________
#####load and parse files
g1.p <- csv_parse(g1.pos)
g6.p <- csv_parse(g6.pos)

g1.n <- csv_parse(g1.neg)
g6.n <- csv_parse(g6.neg)

#Combine HILIC pos and neg groups together
full.pos <- bind_rows(g1.p, g6.p) %>%
  mutate_if(.,is.character,str_replace_all, pattern="SBf", replacement = "SBe")
full.neg <- bind_rows(g1.n, g6.n) %>%
  mutate_if(.,is.character,str_replace_all, pattern="SBf", replacement = "SBe")

###remove standard mixes and blanks 
# from the data frame; correct problematic sample names; make dates consistent;

full.pos.ns <- full.pos %>%
  filter(!str_detect(.$Rep, "Std_"))  %>%
  filter(!str_detect(.$Rep, "LinA"))  %>%
  filter(!str_detect(.$Rep, "blk")) %>%
  mutate(Rep = str_replace_all(Rep, pattern = "Smp_A_SBe_Aa", replacement = "Smp_A_SBe_A")) %>%
  mutate(Rep = str_replace_all(Rep, pattern = "200914_Smp_PS_NS_A", replacement = "Smp_PS_NS_A")) %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point1_S_", replacement = "LinA_0point1_") %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point5_S_", replacement = "LinA_0point5_") %>%
  mutate(date = "200914_Smp") %>%
  unite("Rep", c("date", "Rep")) %>%
  mutate_if(.,is.character,str_replace_all, pattern="Smp_Poo", replacement = "Poo")

full.neg.ns <- full.neg %>%
  filter(!str_detect(.$Rep, "Std_"))  %>%
  filter(!str_detect(.$Rep, "LinA"))  %>%
  filter(!str_detect(.$Rep, "blk")) %>%
  mutate(Rep = str_replace_all(Rep, pattern = "Smp_A_SBe_Aa", replacement = "Smp_A_SBe_A")) %>%
  mutate(Rep = str_replace_all(Rep, pattern = "200914_Smp_PS_NS_A", replacement = "Smp_PS_NS_A")) %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point1_S_", replacement = "LinA_0point1_") %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point5_S_", replacement = "LinA_0point5_") %>%
  mutate(date = "200914_Smp") %>%
  unite("Rep", c("date", "Rep"))  %>%
  mutate_if(.,is.character,str_replace_all, pattern="Smp_Poo", replacement = "Poo")


#####Make and export IS files
pos.IS <- full.pos.ns %>%
  filter(str_detect(.$Compound, ","))
write_csv(pos.IS, file = "Intermediates/Analytical_Validation/AV_PPL_HILIC_Pos_IS.csv")

neg.IS <- full.neg.ns %>%
  filter(str_detect(.$Compound, ","))
write_csv(neg.IS, file = "Intermediates/Analytical_Validation/AV_PPL_HILIC_Neg_IS.csv")


####Make and export data file for BMIS Normalization by removing internal standards, 
## and adding in a polarity identifier to compound names and then combining pos and neg 
# data into a single data frame and then pivot wider to prepare for BMIS code
pos.AV <- full.pos.ns %>%
  filter(!str_detect(.$Compound, ",")) %>%
  mutate(Mode = "Pos") %>%
  unite("Compound", c("Compound", "Mode"))

neg.AV <- full.neg.ns %>%
  filter(!str_detect(.$Compound, ",")) %>%
  mutate(Mode = "Neg") %>%
  unite("Compound", c("Compound", "Mode"))

full.AV <- bind_rows(pos.AV, neg.AV)
full.AV.export <- pivot_wider(full.AV, names_from = Rep, values_from = Area)
write_csv(full.AV.export, file = "Intermediates/Analytical_Validation/AV_PPL_HILIC_combined_raw.csv")


#######Clean up Sample List associated with data set and write to a new file
AV.smp.list <- read_csv(smp.list, skip = 1) %>%
  rename(Rep = `File Name`) %>%
  filter(str_detect(.$Rep, "PPL") |
           str_detect(.$Rep, "Poo")) %>%
  filter(!str_detect(.$Rep, "Blk")) %>%
  filter(!str_detect(.$Rep, "blk")) %>%
  mutate_if(.,is.character,str_replace_all, pattern="PPL_", replacement = "") %>%
  mutate_if(.,is.character,str_replace_all, pattern="200914_", replacement = "") %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point1_S_", replacement = "LinA_0point1_") %>%
  mutate_if(.,is.character,str_replace_all, pattern="LinA_0point5_S_", replacement = "LinA_0point5_") %>%
  filter(!str_detect(.$Rep, "NoSpike"))%>%
  mutate(date = "200914") %>%
  unite("Rep", c("date", "Rep")) %>%
  rename("Injec_vol" = InjVol) %>%
  mutate_if(.,is.character,str_replace_all, pattern="SBf", replacement = "SBe")

write_csv(AV.smp.list, file = "Intermediates/Analytical_Validation/AV_PPL_HILIC_Samp_List.csv")

