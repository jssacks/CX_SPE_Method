

#library
library(readr)
library(tidyverse)
####Quality Control Script

##Define inputs:
smp.file <- "Intermediates/Environmental_Samples/ES_CX_RP_BMISed_dat.csv"
blk.file <- "Intermediates/Environmental_Samples/ES_CX_RP_targeted_QCflags.csv"

#####Blanks QC
Sample.dat <- read_csv(smp.file) 
Blank.dat <- read_csv(blk.file) %>%
  filter(str_detect(.$Rep, "Blk")) %>%
  filter(!str_detect(.$Rep, "PPL")) %>%
  filter(!str_detect(.$Rep, "Blank6")) %>%    
  filter(!str_detect(.$Rep, "Blank7")) %>%    
  mutate(Area = replace_na(Area, 0)) 

Blk.ave.dat <- Blank.dat %>%
  group_by(Compound) %>%
  summarize(Blk.Av = mean(Area),
            Blk.sd = sd(Area),
            Blk.max = max(Area),
            Blk.LD = Blk.Av + (3.182 * (Blk.sd/sqrt(14)))) %>%
  rename(MF = Compound)

write_csv(Blk.ave.dat, file = "Intermediates/Environmental_Samples/ES_CX_RP_Blk_LOD_signal_values.csv")

Blk.QC <- left_join(Sample.dat, Blk.ave.dat) %>%
  mutate(Blk.flag = ifelse(Adjusted_Area <= Blk.LD, 1, 0))
Blk.QC.list <- Blk.QC %>%
  group_by(MF, samp) %>%
  summarise(Blk.flag = sum(Blk.flag)) %>%
  filter(!Blk.flag > 0) %>%
  mutate(Blk.QC = "ok") %>%
  select(MF, samp, Blk.QC)

Sample.dat.blkQC <- left_join(Sample.dat, Blk.QC.list) %>%
  filter(Blk.QC == "ok")

####QC flags 
flag.dat <- read_csv(blk.file) %>%
  rename(MF = Compound)
sample.flag.dat <- left_join(Sample.dat.blkQC, flag.dat)

sample.flag.dat.remove <- sample.flag.dat %>%
  filter(is.na(areaminFlag))

write_csv(sample.flag.dat.remove, file = "Intermediates/Environmental_Samples/ES_CX_RP_BMISed_BlkQCed_dat.csv")
