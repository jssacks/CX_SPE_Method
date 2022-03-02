#
#
#

#
##
#
#
#
#
library(tidyverse)
library(readr)
####
#
#
#define inputs
CX.HQcomps.file <- "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv"
ANOVA.file <- "Intermediates/Analytical_Validation/AV_CX_ANOVA_Results.csv"

###
CX.dat <- read_csv(CX.HQcomps.file)
ANOVA.dat <- read_csv(ANOVA.file)
