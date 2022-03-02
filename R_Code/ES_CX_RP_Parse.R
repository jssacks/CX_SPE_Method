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
g1.pos <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_RP_G1.csv"

#sample list
smp.list <- "Meta_Data/Sample_Lists/CXSPE_ES_RP_SampList.csv"


#################### ___________________________________________
#####load and parse files
full.pos <- ES_csv_parse(g1.pos)

###QC
###Quality Control Flags for Area, SN, and Mass Error
###QC Parameters:
SNmin = 1
ppmflex = 6
Areamin = 5000

###Pos QC
full.pos.qc <- full.pos %>%
  mutate(SN = Area/Background) %>%
  mutate(SNFlag = ifelse((Area/Background < SNmin), "SNFlag", NA)) %>%
  mutate(ppmFlag = ifelse((abs(Mass.Error.PPM) > ppmflex), "ppmFlag", NA)) %>%
  mutate(areaminFlag = ifelse((Area < Areamin), "areaminFlag", NA))

###Write QC Flag data to CSV
qc.export <- full.pos.qc %>%
  filter(!str_detect(.$Compound, ", "))
write_csv(qc.export, path = "Intermediates/Environmental_Samples/ES_CX_RP_targeted_QCFlags.csv")


####Get Targeted Data and export:
rp.targeted <- full.pos %>%
  filter(!str_detect(.$Compound, ", "))
rp.targeted.export <- pivot_wider(id_cols = Compound, rp.targeted, names_from = Rep, values_from = Area)
write_csv(rp.targeted.export, path = "Intermediates/Environmental_Samples/ES_CX_RP_targeted_combined_raw.csv")

###Get IS Data and export:
rp.IS <- full.pos %>%
  filter(str_detect(.$Compound, ", "))
write_csv(rp.IS, path = "Intermediates/Environmental_Samples/ES_CX_RP_IS.csv")


###Figure out Sample Lists
samp.list <- read_csv(smp.list,
                      skip = 1)
halfpoos <- samp.list %>%
  filter(str_detect(.$`File Name`, "Half")) %>%
  mutate(InjVol = InjVol/2)
nonpoos <- samp.list %>%
  filter(!str_detect(.$`File Name`, "Half"))
samp.list.exp <- rbind(halfpoos, nonpoos) %>%
  mutate(Rep = `File Name`) %>%
  select(Rep, InjVol)
write_csv(samp.list.exp, path = "Intermediates/Environmental_Samples/ES_RP_SampList.csv")



































