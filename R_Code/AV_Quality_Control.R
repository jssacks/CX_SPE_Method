##
#
#
#This script...
#
#

#packages
library(tidyverse)
library(readr)


###Define inputs
CX.EE.file <- "Intermediates/Analytical_Validation/AV_CX_EE_RF_Dat.csv"
CX.R2.file <- "Intermediates/Analytical_Validation/AV_CX_R2_Dat.csv"

PPL.EE.file <- "Intermediates/Analytical_Validation/AV_PPL_EE_RF_Dat.csv"

###Define Quality Control Cutoffs:
EE.cut.high <- 150
EE.cut.low <- 1
R2.cut <- 0.7
RSD.cut <- 50
#RR.cut <-  50

###Apply quality control cutoffs to CX-SPE ALOHA Data:
CX.EE.dat <- read_csv(CX.EE.file)
CX.R2.dat <- read_csv(CX.R2.file)

CX.dat <- left_join(CX.EE.dat, CX.R2.dat) %>%
  filter(sample == "A") %>%
  mutate(Overall.Mean.EE = as.numeric(Overall.Mean.EE),
         Overall.SD.EE = as.numeric(Overall.SD.EE),
         Sample.Mean.EE = as.numeric(Sample.Mean.EE),
         Sample.SD.EE = as.numeric(Sample.SD.EE),
         Sample.RSD = as.numeric((Sample.SD.EE/Sample.Mean.EE)*100),
         Overall.RSD = as.numeric((Overall.SD.EE/Overall.Mean.EE)*100),
         Sample.RR = as.numeric((Sample.Range.EE/Overall.Mean.EE)*100),
         R2 = as.numeric(R2)) %>%
  select(MF, Fraction, sample, Overall.Mean.EE, Overall.RSD, 
         Overall.SD.EE, Sample.Mean.EE, Sample.SD.EE, Sample.Range.EE, Sample.RR,
         Sample.RSD, R2, z, column) %>%
  unique()

CX.A.QC <- CX.dat%>%
  filter(as.numeric(Overall.Mean.EE) <= EE.cut.high) %>%
  filter(as.numeric(Overall.Mean.EE) >= EE.cut.low) %>%
  filter(as.numeric(R2) >= R2.cut) %>%
  filter(as.numeric(Overall.RSD) <= RSD.cut)
#  filter(as.numeric(Sample.RR) <= RR.cut)

####Make leucine and isoleucine into (Iso)leucine 
CX.A.QC.2 <- CX.A.QC %>%
  filter(!MF == "Leucine") %>%
  mutate_if(.,is.character,str_replace_all, pattern="Isoleucine", replacement = "(Iso)leucine")





####
write_csv(CX.A.QC.2, file = "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv")


###Apply quality control cutoffs to PPL ALOHA Data
PPL.EE.dat <- read_csv(PPL.EE.file)

PPL.dat <- PPL.EE.dat %>%
  filter(sample == "A") %>%
  mutate(Overall.Mean.EE = as.numeric(Overall.Mean.EE),
         Overall.SD.EE = as.numeric(Overall.SD.EE),
         Sample.Mean.EE = as.numeric(Sample.Mean.EE),
         Sample.SD.EE = as.numeric(Sample.SD.EE),
         Sample.RSD = as.numeric((Sample.SD.EE/Sample.Mean.EE)*100),
         Overall.RSD = as.numeric((Overall.SD.EE/Overall.Mean.EE)*100)) %>%
  select(MF, Fraction, sample, Overall.Mean.EE, Overall.RSD,
         Overall.SD.EE, Sample.Mean.EE, Sample.SD.EE, 
         Sample.RSD) %>%
  unique()

PPL.A.QC <- PPL.dat%>%
  filter(as.numeric(Overall.Mean.EE) <= EE.cut.high) %>%
  filter(as.numeric(Overall.Mean.EE) >= EE.cut.low) %>%
  filter(as.numeric(Overall.RSD) <= RSD.cut)

###Make leucine and isoleucine into (iso)leucine
PPL.A.QC.2 <- PPL.A.QC %>%
  filter(!MF == "Leucine") %>%
  mutate_if(.,is.character,str_replace_all, pattern="Isoleucine", replacement = "(Iso)leucine")


####
write_csv(PPL.A.QC.2, file = "Intermediates/Analytical_Validation/AV_PPL_HQ_dat.csv")



















