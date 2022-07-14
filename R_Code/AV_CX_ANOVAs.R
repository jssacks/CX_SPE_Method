##
#
#
#This script...
#
#

#packages
library(tidyverse)
library(readr)

#Source Functions
source("R_Code/CXSPE_Functions.R")

###Define inputs
CX.EE.file <- "Intermediates/Analytical_Validation/AV_CX_EE_RF_Dat.csv"

CX.HQcomps.file <- "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv"

###
CX.dat <- read_csv(CX.EE.file) %>%
  filter(!MF == "Leucine") %>%
  mutate_if(.,is.character,str_replace_all, pattern="Isoleucine", replacement = "(Iso)leucine")

CX.HQcomps <- read_csv(CX.HQcomps.file) %>%
  select(MF)

CX.HQ.dat <- left_join(CX.HQcomps, CX.dat)

####Run ANOVA to test for differences in EE across samples for each HQ compound

#create output dataframe
out.EE <- data.frame(
  "MF" = character(),
  "Pval.EE" = numeric())

#create MF list to put through loop
dat.MF <- CX.HQ.dat %>%
  ungroup()%>%
  select(MF) %>%
  unique()

#run loop for EE_anova function for all HQ CX compounds
for (i in seq_along(dat.MF$MF)) {
  res <- dat.MF$MF[i]
  EE.stat <- EE_anova(CX.HQ.dat, res)
  res.df <- data.frame(EE.stat)
  out.EE <- rbind(out.EE, res.df) %>%
    unique()
}

#Apply FDR p-value adjustment 
out.EE$Pval.EE <- p.adjust(out.EE$Pval.EE, method = "BH")

write_csv(out.EE, file = "Intermediates/")

###Identify compounds that have significant differences in EE (FDR adjusted p-value <0.05)
EE.sig <- full_join(CX.HQ.dat, out.EE) %>%
  mutate(Signif.EE = ifelse(Pval.EE <(0.05), "Yes", "No"))




####Run ANOVA to test for differences in RF across samples for each HQ compound

#create output dataframe
out.RF <- data.frame(
  "MF" = character(),
  "Pval.rf" = numeric())

#####
for (i in seq_along(dat.MF$MF)) {
  res <- dat.MF$MF[i]
  EE.stat <- RF_anova(CX.HQ.dat, res)
  res.df <- data.frame(EE.stat)
  out.RF <- rbind(out.RF, res.df) %>%
    unique()
}

#Apply FDR p-value adjustment 
out.RF$Pval.rf <- p.adjust(out.RF$Pval.rf, method = "BH")

###Identify compounds that have significant differences in EE (FDR adjusted p-value <0.05)
EE.RF.sig <- full_join(EE.sig, out.RF) %>%
  mutate(Signif.RF = ifelse(Pval.rf <(0.05), "Yes", "No"))


###write EE and RF output to a csv file
write_csv(EE.RF.sig, file = "Intermediates/Analytical_Validation/AV_CX_ANOVA_Results.csv")





















