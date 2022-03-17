#
#
#
#
#
#This script...

#load packages
library(tidyverse)
library(readr)
library(stringr)

#load functions
source("R_Code/CXSPE_Functions.R")


#Define inputs
#HILIC positive skyline ES files for CX-SPE
g1.pos <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_Pos_G1.csv"
g2.pos <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_Pos_G2.csv"
g3.pos <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_Pos_G3.csv"

#HILIC negative skyline ES files for CX-SPE
g1.neg <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_Neg_G1.csv"
g2.neg <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_Neg_G2.csv"
g3.neg <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_Neg_G3.csv"

#sample list
smp.list <- "Meta_Data/Sample_Lists/CXSPE_ES_HILIC_SampList.csv"


#################### ___________________________________________
#####load and parse files
g1.p <- ES_csv_parse(g1.pos)
g2.p <- ES_csv_parse(g2.pos)
g3.p <- ES_csv_parse(g3.pos)
#
g1.n <- ES_csv_parse(g1.neg)
g2.n <- ES_csv_parse(g2.neg)
g3.n <- ES_csv_parse(g3.neg)


#Combine all HILIC groups together
full.pos <- bind_rows(g1.p, g2.p) %>%
  filter(!str_detect(.$Rep, "Std")) %>%
  filter(!str_detect(.$Rep, "Mort")) %>%
  bind_rows(., g3.p)
full.neg <- bind_rows(g1.n, g2.n)  %>%
  filter(!str_detect(.$Rep, "Std"))  %>%
  filter(!str_detect(.$Rep, "Mort")) %>%
  bind_rows(., g3.n)

###Quality Control Flags for Area, SN, and Mass Error
###QC Parameters:
SNmin = 1
ppmflex = 6
Areamin = 40000

###Pos QC
full.pos.qc <- full.pos %>%
  mutate(SN = Area/Background) %>%
  mutate(SNFlag = ifelse((Area/Background < SNmin), "SNFlag", NA)) %>%
  mutate(ppmFlag = ifelse((abs(Mass.Error.PPM) > ppmflex), "ppmFlag", NA)) %>%
  mutate(areaminFlag = ifelse((Area < Areamin), "areaminFlag", NA))

full.neg.qc <- full.neg %>%
  mutate(SNFlag = ifelse((Area/Background < SNmin), "SNFlag", NA)) %>%
  mutate(ppmFlag = ifelse((abs(Mass.Error.PPM) > ppmflex), "ppmFlag", NA)) %>%
  mutate(areaminFlag = ifelse((Area < Areamin), "areaminFlag", NA))




###add in ion mode and date for targeted data
full.pos.targeted <- full.pos.qc %>%
  filter(!str_detect(.$Compound, ",")) %>%
  mutate(Mode = "Pos") %>%
  unite("Compound", c("Compound", "Mode")) %>%
  mutate(date = "210322") %>%
  unite("Rep", c("date", "Rep"))

full.neg.targeted <- full.neg.qc %>%
  filter(!str_detect(.$Compound, ",")) %>%
  mutate(Mode = "Neg") %>%
  unite("Compound", c("Compound", "Mode")) %>%
  mutate(date = "210322") %>%
  unite("Rep", c("date", "Rep"))

###add in ion mode and date for IS data
ES.full.pos.IS <- full.pos %>%
  filter(str_detect(.$Compound, ",")) %>%
  mutate(date = "210322") %>%
  unite("Rep", c("date", "Rep")) %>%
  select(Rep, Compound, Area)
write_csv(ES.full.pos.IS, file = "Intermediates/Environmental_Samples/ES_HILIC_pos_IS.csv")


ES.full.neg.IS <- full.neg %>%
  filter(str_detect(.$Compound, ",")) %>%
  mutate(date = "210322") %>%
  unite("Rep", c("date", "Rep")) %>%
  select(Rep, Compound, Area)
write_csv(ES.full.neg.IS, file = "Intermediates/Environmental_Samples/ES_HILIC_neg_IS.csv")


###Combine targeted data and export for BMIS and for later QC
ES.full.HILIC.targeted <- bind_rows(full.pos.targeted, full.neg.targeted)

###write QC flag data to csv:
write_csv(ES.full.HILIC.targeted, file = "Intermediates/Environmental_Samples/ES_HILIC_targeted_combined_QCflags.csv")

###write data for BMIS
ES.full.HILIC.targeted.export <- ES.full.HILIC.targeted %>%
  select(Rep, Compound, Area)
ES.full.HILIC.targeted.export.2 <- pivot_wider(ES.full.HILIC.targeted.export, names_from = Rep, values_from = Area)
write_csv(ES.full.HILIC.targeted.export.2, file = "Intermediates/Environmental_Samples/ES_HILIC_targeted_combined_raw.csv")


###Figure out Sample Lists
samp.list <- read_csv(smp.list, skip = 1)
halfpoos <- samp.list %>%
  filter(str_detect(.$`File Name`, "Half")) %>%
  mutate(InjVol = InjVol/2)
nonpoos <- samp.list %>%
  filter(!str_detect(.$`File Name`, "Half"))
samp.list.exp <- rbind(halfpoos, nonpoos)
write_csv(samp.list.exp, file = "Intermediates/Environmental_Samples/ES_HILIC_SampList.csv")

















