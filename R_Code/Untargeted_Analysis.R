#####
#
#
#
#
#
#this script....
#
#
#

#
library(tidyverse)
library(readr)
library(ggpubr)

###input
long.dat.file <- "Intermediates/Untargeted/Untargeted_QCed_dat_long.csv"

###
dat <- read_csv(long.dat.file)

####get just 1 value for each triplicate
dat.1 <- dat %>%
  select(MF, Fraction, method, Fraction, Retention_time, mz, Station) %>%
  unique()


#####Summarize number of MFs in each fraction
dat.sum.frac <- dat.1 %>%
  select(MF, Fraction) %>%
  unique() %>%
  group_by(Fraction) %>%
  summarise(MFs_observed = n())
dat.sum.frac

#####Summarize number of MFs in each sample
dat.sum.sta <- dat.1 %>%
  select(MF, Station) %>%
  unique() %>%
  group_by(Station) %>%
  summarise(MFs_observed = n())
dat.sum.sta


###Statistics

##tidy up data
dat.2 <- dat.1 %>%
  mutate(chrom.col = case_when(str_detect(.$Fraction, "HILIC") ~ "HILIC",
                               str_detect(.$Fraction, "RP") ~ "RP")) %>%
  mutate(method = str_replace(.$method, "CatEx", "CX-SPE")) %>%
  mutate(method = str_replace(.$method, "PPL", "PPL-SPE")) %>%
  mutate(Station = str_replace_all(.$Station, "A", "ALOHA"))


#####HILIC Comparisons_________________________________________________________
dat.2.HILIC <- dat.2 %>%
  filter(chrom.col == "HILIC") 
dat.2.HILIC.CX <- dat.2.HILIC %>%
  filter(method == "CX-SPE")
dat.2.HILIC.PPL <- dat.2.HILIC %>%
  filter(method == "PPL-SPE")

#####Check to see if RT of samples (ALOHA and PS) are significantly different 
kruskal.test(Retention_time ~ Station, dat.2.HILIC.CX)
kruskal.test(Retention_time ~ Station, dat.2.HILIC.PPL)

#####Check to see if m/z of samples (ALOHA and PS) are significantly different 
kruskal.test(mz ~ Station, dat.2.HILIC.CX)
kruskal.test(mz ~ Station, dat.2.HILIC.PPL)

######Compare CX-SPE and PPL-SPE by RT
kruskal.test(Retention_time ~ method, dat.2.HILIC)

######Compare CX-SPE and PPL-SPE by m/z
kruskal.test(mz ~ method, dat.2.HILIC)



#####RP Comparisons_____________________________________________________________
dat.2.RP <- dat.2 %>%
  filter(chrom.col == "RP") 
dat.2.RP.CX <- dat.2.RP %>%
  filter(method == "CX-SPE")
dat.2.RP.PPL <- dat.2.RP %>%
  filter(method == "PPL-SPE")
dat.2.RP.A <- dat.2.RP %>%
  filter(Station == "ALOHA")
dat.2.RP.PS <- dat.2.RP %>%
  filter(Station == "PS")

#See if stations are different with respect to RT
kruskal.test(Retention_time ~ Station, dat.2.RP.CX)
kruskal.test(Retention_time ~ Station, dat.2.RP.PPL)

#_______since stations are different, Compare RT for ALOHA and PS separately 
#Aloha
kruskal.test(Retention_time ~ method, dat.2.RP.A)
#Puget Sound
kruskal.test(Retention_time ~ method, dat.2.RP.PS)


###See if stations are different with respect to m/z 
kruskal.test(mz ~ Station, dat.2.RP.CX)
kruskal.test(mz ~ Station, dat.2.RP.PPL)


#_______since stations are different, Compare m/z for ALOHA and PS separately 
#Aloha
kruskal.test(mz ~ method, dat.2.RP.A)
#Puget Sound
kruskal.test(mz ~ method, dat.2.RP.PS)






