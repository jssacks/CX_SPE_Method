

#this script...



#library
library(readr)
library(tidyverse)
####Quality Control Script

##Define inputs:
smp.file <- "Intermediates/Environmental_Samples/ES_CatEx_HILIC_BMISed_dat.csv"
blk.file <- "Intermediates/Environmental_Samples/ES_HILIC_targeted_combined_QCflags.csv"

###define calculation inputs:

#number of blanks (n)
n.blks <- 15

#students t-value for n-1 at 95% confidence interval
tval <- 1.761




#####Blanks QC
Sample.dat <- read_csv(smp.file) 
Blank.dat <- read_csv(blk.file) %>%
  filter(str_detect(.$Rep, "Blk")) %>%
  filter(!str_detect(.$Rep, "PPL")) %>%
  filter(!str_detect(.$Rep, "Blank6")) %>%    ####NEED TO CHECK ON THIS
  filter(!str_detect(.$Rep, "Blank7")) %>%    
  mutate(Area = replace_na(Area, 0)) 

Blk.ave.dat <- Blank.dat %>%
  group_by(Compound) %>%
  summarize(Blk.Av = mean(Area),
            Blk.sd = sd(Area),
            Blk.max = max(Area),
            Blk.LD = Blk.Av + (tval * (Blk.sd/sqrt(n.blks)))) %>%
  rename(MF = Compound)

####Address issue with Cytidine, Ribolflavin Monophosphate, and Inosine having 0 signal in blanks 
# by setting blank LD value at 1/2 minimum area threshold for a peak (20,000)
Blk.ave.dat.corrected <- Blk.ave.dat %>%
  mutate(Blk.LD = case_when(
    .$MF == "Cytidine_Pos" ~ 20000,
    .$MF == "Riboflavin Monophosphate_Pos" ~ 20000,
    .$MF == "Inosine_Neg" ~ 20000,
    TRUE ~ Blk.LD
  ))


write_csv(Blk.ave.dat.corrected, file = "Intermediates/Environmental_Samples/ES_CX_HILIC_Blk_LOD_signal_values.csv")

Blk.QC <- left_join(Sample.dat, Blk.ave.dat.corrected) %>%
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
flag.dat <- read_csv("Intermediates/Environmental_Samples/ES_HILIC_targeted_combined_QCflags.csv") %>%
  rename(MF = Compound)
sample.flag.dat <- left_join(Sample.dat.blkQC, flag.dat)

sample.flag.dat.remove <- sample.flag.dat %>%
  filter(is.na(areaminFlag))

write_csv(sample.flag.dat.remove, file = "Intermediates/Environmental_Samples/ES_CatEx_HILIC_BMISed_BlkQCed_dat.csv")
