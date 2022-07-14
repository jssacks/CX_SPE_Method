
####Quality Control Script

#####Blanks QC
Sample.dat <- read_csv("Intermediates/Environmental_Samples/ES_PPL_HILIC_BMISed_dat.csv") 
#Blank.dat <- read_csv("Intermediates/Environmental_Samples/ES_HILIC_targeted_combined_QCflags.csv") %>%
#  pivot_longer(cols = -`Compound`, names_to = "Rep", values_to = "Area") %>%

#number of blanks (n)
n.blks <- 3

#students t-value for n-1 at 95% confidence interval
tval <- 2.353


Blank.dat <- Sample.dat %>%
  filter(str_detect(.$SampID, "Blk")) %>%
  mutate(Area = replace_na(Area, 0)) %>%
  rename(Compound = MF)


Blk.ave.dat <- Blank.dat %>%
  group_by(Compound) %>%
  summarize(Blk.Av = mean(Area),
            Blk.sd = sd(Area),
            Blk.max = max(Area),
            Blk.LD = Blk.Av + (tval * (Blk.sd/sqrt(n.blks)))) %>%
  rename(MF = Compound)

Blk.ave.dat.corrected <- Blk.ave.dat %>%
  mutate(Blk.LD = case_when(
    .$Blk.LD == 0 ~ 20000,
    is.na(.$Blk.LD) ~ 20000, 
    TRUE ~ Blk.LD
  ))


write_csv(Blk.ave.dat.corrected, file = "Intermediates/Environmental_Samples/ES_PPL_Blk_LOD_signal_values.csv")
Sample.dat.2 <- Sample.dat %>%
  filter(!str_detect(.$SampID, "Blk")) %>%
  filter(!str_detect(.$SampID, "Poo"))


Blk.QC <- left_join(Sample.dat.2, Blk.ave.dat) %>%
  mutate(Blk.flag = ifelse(Adjusted_Area <= Blk.LD, 1, 0))
Blk.QC.list <- Blk.QC %>%
  group_by(MF, samp) %>%
  summarise(Blk.flag = sum(Blk.flag)) %>%
  filter(!Blk.flag > 0) %>%
  mutate(Blk.QC = "ok") %>%
  select(MF, samp, Blk.QC)

Sample.dat.blkQC <- left_join(Sample.dat.2, Blk.QC.list) %>%
  filter(Blk.QC == "ok")

####QC flags 
flag.dat <- read_csv("Intermediates/Environmental_Samples/ES_HILIC_targeted_combined_QCflags.csv") %>%
  rename(MF = Compound)
sample.flag.dat <- left_join(Sample.dat.blkQC, flag.dat)

sample.flag.dat.remove <- sample.flag.dat %>%
  filter(is.na(areaminFlag))

write_csv(sample.flag.dat.remove, file = "Intermediates/Environmental_Samples/ES_PPL_HILIC_BMISed_BlkQCed_dat.csv")

