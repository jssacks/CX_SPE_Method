library(ggplot2)
library(tidyverse)
options(dplyr.summarise.inform=F)
options(readr.num_columns = 0)

#Define all your inputs here------
HILIC.longdat.file <- "Intermediates/Untargeted/Untargeted_HILIC_BMIS_Dat.csv"
CyanoAq.longdat.file  <- "Intermediates/Untargeted/Untargeted_RP_BMIS_Dat.csv"
smp.to.blank.ratio <- 3
rsd.of.poo.max <- 0.3
min.peak.area <- 5000
min.SN <- 5
HILIC.SN.file.1 <- "Raw_Data/MS_Dial_Output/HILIC_Pos/SN_2_2021818153.txt"
HILIC.SN.file.2 <- "Raw_Data/MS_Dial_Output/HILIC_Neg/SN_0_20218181538.txt"
CyanoAq.SN.file <- "Raw_Data/MS_Dial_Output/RP/SN_0_2021818167.txt"
#Vol.filt <- "Sample_Lists/MEP_volume_filtered.csv"



#Load and clean up the SNs-----
SN.dat.HILIC.1 <- read_delim(HILIC.SN.file.1,
                             "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4) %>%
  mutate(Column = "HILICPos") %>% 
  mutate(MF = paste0(`Alignment ID`, "_", Column)) %>%
  rename(Retention_time = `Average Rt(min)`,
         mz = `Average Mz`,
         annotation = `Metabolite name`) %>%
  select(MF, Retention_time, mz, annotation, `210322_Blk_Blank1_A`:`210322_Smp_PS_C`) %>%
  pivot_longer(`210322_Blk_Blank1_A`:`210322_Smp_PS_C`, names_to = "SampID", values_to = "SN") 

SN.dat.HILIC.2 <- read_delim(HILIC.SN.file.2,
                             "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4) %>%
  mutate(Column = "HILICNeg") %>% 
  mutate(MF = paste0(`Alignment ID`, "_", Column)) %>%
  rename(Retention_time = `Average Rt(min)`,
         mz = `Average Mz`,
         annotation = `Metabolite name`) %>%
  select(MF, Retention_time, mz, annotation, `210322_Blk_Blank1_A`:`210322_Smp_PS_C`) %>%
  pivot_longer(`210322_Blk_Blank1_A`:`210322_Smp_PS_C`, names_to = "SampID", values_to = "SN")

SN.dat.CyanoAq <- read_delim(CyanoAq.SN.file,
                             "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4) %>%
  mutate(Column = "RP") %>% 
  mutate(MF = paste0(`Alignment ID`, "_", Column)) %>%
  rename(Retention_time = `Average Rt(min)`,
         mz = `Average Mz`,
         annotation = `Metabolite name`) %>%
  select(MF, Retention_time, mz, annotation, `210329_Blk_Blank1_A`:`210329_Smp_PS_C`) %>%
  pivot_longer(`210329_Blk_Blank1_A`:`210329_Smp_PS_C`, names_to = "SampID", values_to = "SN") 



sn.dat.combo <-SN.dat.HILIC.1 %>% bind_rows(SN.dat.HILIC.2) %>%
  bind_rows(SN.dat.CyanoAq) 

MF.info.combo <-sn.dat.combo %>%
  select(MF, Retention_time, mz, annotation) %>% unique()


#Load your BMISd longdats match by samples, attach the SN data------
HILIC.BMISd.dat <- read_csv(HILIC.longdat.file)
RP.BMISd.dat <- read_csv(CyanoAq.longdat.file)

dat <- HILIC.BMISd.dat %>%
  bind_rows(RP.BMISd.dat)%>%
  left_join(sn.dat.combo, by = c("MF", "SampID")) %>%
  filter(!str_detect(.$SampID, "DDA")) %>%
  filter(!str_detect(.$SampID, "Filter_blk_2")) %>%
  mutate(method = case_when(str_detect(.$SampID, "PPL") ~ "PPL",
                            !str_detect(.$SampID, "PPL") ~ "CatEx")) %>%
  mutate(Station = case_when(str_detect(.$SampID, "_Aloha_") ~ "A",
                             str_detect(.$SampID, "_PS_") ~ "PS"))

#For each MF, get the following information: sample type blank area; CV of pooled (per sample type)-----
dat.RSD.poo <- dat %>%
  filter(type == "Poo") %>%
  group_by(MF) %>%
  summarise(RSD_poo = sd(Adjusted_Area)/mean(Adjusted_Area)) 

dat.blank.ave <- dat %>%
  filter(type == "Blk") %>%
  group_by(MF, method) %>%
  summarise(ave_blank = mean(Adjusted_Area)) %>%
  mutate(ave_blank = ifelse(ave_blank == 0, 100, ave_blank))

dat.smp.ave <- dat %>%
  filter(type == "Smp") %>%
  group_by(MF, method, Station) %>%
  summarise(ave_smp = mean(Adjusted_Area)) 

dat.MF.good <- dat.RSD.poo %>%
  left_join(dat.blank.ave, by = c("MF")) %>%
  left_join(dat.smp.ave, by = c("MF", "method")) %>%
  mutate(smptoblank = ave_smp/ave_blank) %>%
  mutate(GoodMFs = ifelse(smptoblank > smp.to.blank.ratio & RSD_poo < 0.3, 1, 0)) %>%
  mutate(GoodMFs = ifelse(is.na(GoodMFs), 0, GoodMFs)) 

#See if each individual peaks are good enough by minimum peak size and maximum signal to noise-----
dat.QC <- dat %>%
  mutate(GoodPeak = ifelse(Adjusted_Area > min.peak.area & SN > 5, 1, 0)) %>%
  left_join(dat.MF.good, by = c("MF", "method", "Station"))


#Summarize number of "good mass features in each fraction and in each sample type"-----
dat.QC.summary <- dat.QC %>%
  mutate(Fraction = str_extract(MF, "_\\w+$") %>%
           str_replace("_", "")) %>%
  filter(type == "Smp") %>%
  select(MF, SampID:GoodPeak, GoodMFs, Fraction, method, Station) 
dat.QC.summary.2 <- dat.QC.summary %>%
  group_by(MF, Fraction, method, Station) %>%
  summarise(GoodPeakTotal = sum(GoodPeak),
            GoodMFsTotal = sum(GoodMFs)) %>%
  mutate(Keep = ifelse(GoodPeakTotal >= 3 & GoodMFsTotal >= 3, 1, 0)) 
test <- dat.QC.summary.2 %>% filter(is.na(Keep))

dat.QC.summary.3 <- dat.QC.summary.2 %>%
  group_by(Fraction, method, Station) %>%
  summarise(MFs_observed = sum(Keep))
#write_csv(dat.QC.summary.3, "MF_observed_summary.csv")


####
dat.keep.id <- dat.QC.summary.2 %>%
  filter(Keep == 1) %>%
  select(MF, Fraction, method, Station, Keep)



#Make supplemental table
dat.QC.2 <- left_join(dat.keep.id, dat.QC) %>%
  rename(QCd_area = Adjusted_Area) %>%
  select(MF, SampID, Retention_time, mz, QCd_area, method, Station)

full.dat.QC <- dat.QC.2

#make long dat file
write_csv(full.dat.QC, path = "Intermediates/Untargeted/Untargeted_QCed_dat_long.csv")


#make wide dat file
dat.QC.2.wide <- full.dat.QC %>%
  select(MF, SampID, QCd_area) %>%
  pivot_wider(names_from = SampID, values_from = QCd_area) %>%
  left_join(MF.info.combo) %>%
  mutate(Fraction = str_extract(MF, "_\\w+$") %>%
           str_replace("_", "")) %>%
  select(MF, Fraction, mz, Retention_time, annotation, everything()) %>%
  arrange(Fraction, mz)

write_csv(dat.QC.2.wide, path = "Intermediates/Untargeted/Untargeted_QCed_dat_wide.csv")
