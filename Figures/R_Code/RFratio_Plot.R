#
#
#
#
#
library(tidyverse)
library(readr)
##
#define inputs
RFratio.file.HILIC <- "Intermediates/Environmental_Samples/HILIC_RFs_RFratios.csv"
RFratio.file.RP <- "Intermediates/Environmental_Samples/ES_CX_RP_RFs_RFratios.csv"
CX.HQ.comps.file <- 
##
HILIC.RF.dat <- read_csv(RFratio.file.HILIC)
RP.RF.dat <- read_csv(RFratio.file.HILIC)