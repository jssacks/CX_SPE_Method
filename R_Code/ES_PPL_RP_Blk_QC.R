####Quality Control Script

#####Blanks QC
Sample.dat <- read_csv("Intermediates/Environmental_Samples/ES_PPL_RP_BMISed_dat.csv") 

Blank.dat <- Sample.dat %>%
  filter(str_detect(.$SampID, "Blk")) %>%
  mutate(Area = replace_na(Area, 0)) %>%
  rename(Compound = MF)

Blk.ave.dat <- Blank.dat %>%
  group_by(Compound) %>%
  summarize(Blk.Av = mean(Area),
            Blk.sd = sd(Area),
            Blk.max = max(Area),
            Blk.LD = Blk.Av + (3.182 * (Blk.sd/sqrt(3)))) %>%
  rename(MF = Compound)

write_csv(Blk.ave.dat, file = "Intermediates/Environmental_Samples/ES_PPL_RP_Blk_LOD_signal_values.csv")

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
flag.dat <- read_csv("Intermediates/Environmental_Samples/ES_CX_RP_targeted_QCflags.csv") %>%
  rename(MF = Compound)
sample.flag.dat <- left_join(Sample.dat.blkQC, flag.dat)

sample.flag.dat.remove <- sample.flag.dat %>%
  filter(is.na(areaminFlag))

write_csv(Sample.dat.blkQC, file  = "Intermediates/Environmental_Samples/ES_PPL_RP_BMISed_BlkQCed_dat.csv")
