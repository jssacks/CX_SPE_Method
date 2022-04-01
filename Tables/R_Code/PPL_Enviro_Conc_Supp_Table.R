#THis script...
#
#
#
library(tidyverse)
library(readr)
###
#
#define input
dat.EE.file <- "Intermediates/Analytical_Validation/AV_PPL_HQ_dat.csv"
dat.HILIC.LOD.file <- "Intermediates/Environmental_Samples/ES_CX_Blk_LOD_Concentrations.csv"
dat.RP.LOD.file <- "Intermediates/Environmental_Samples/ES_CXC_RP_Blk_LOD_Concentrations.csv"
Std.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"

##read in EE data
dat.EE <- read_csv(dat.EE.file) %>%
  rename("Compound" = MF)
