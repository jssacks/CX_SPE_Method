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
CX.HQ.comps.file <- "Intermediates/Analytical_Validation/AV_CX_PPL_HQ_combined.csv"
##


####Load in HILIC and RP RFs and RFratios and combine into a single dataset
HILIC.RF.dat <- read_csv(RFratio.file.HILIC)%>%
  mutate(Compound = str_replace_all(.$Compound, "Isoleucine", "(Iso)leucine")) %>%
  filter(!Compound == "leucine")
#

#Load in RF dat and remove some duplicates
RP.RF.dat <- read_csv(RFratio.file.RP) %>%
  filter(!Compound == "Tyrosine") %>%
  filter(!Compound == "Betaine") %>%
  filter(!Compound == "Butyryl-L-carnitine")

all.RF <- rbind(HILIC.RF.dat, RP.RF.dat) 

##Load in CX-SPE high quality compound details 
CX.HQ.dat <- read_csv(CX.HQ.comps.file) %>%
  select(MF, method, Fraction, sample, Overall.Mean.EE, Overall.SD.EE, 
         Sample.Mean.EE, Sample.SD.EE) %>%
  filter(method == "CX-SPE") %>%
  rename("Compound" = MF) %>%
  mutate(Compound = str_replace_all(.$Compound, "Isoleucine", "(Iso)leucine")) %>%
  filter(!Compound == "leucine")

###Get Just RFratios from CX-SPE HQ Comps
HQ.RF <- left_join(CX.HQ.dat, all.RF) %>%
  unique() %>%
  mutate(RFratio = as.numeric(RFratio)) %>%
  mutate(Fraction = str_replace_all(.$Fraction, "Neg", "HILIC Negative")) %>%
  mutate(Fraction = str_replace_all(.$Fraction, "Pos", "HILIC Positive")) 

###RFratio Plot
rf.plot.1 <- ggplot(HQ.RF, aes(x = RFratio, fill = Fraction)) +
  geom_histogram(color = "black", bins = 40) +
  theme_test() +
  scale_y_continuous(expand = c(0, NA), limits = c(0,25)) +
  xlim(0.25,1.75) +
  geom_vline(xintercept = 1, alpha = 0.7, size = 1.5) +
  geom_vline(xintercept = 1.1, alpha = 0.7, size = 1.5, linetype = "dashed") +
  geom_vline(xintercept = 0.9, alpha = 0.7, size = 1.5, linetype = "dashed") 
rf.plot.1
