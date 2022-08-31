#
#
library(tidyverse)
library(readr)
###


###input
long.dat.file <- "Intermediates/Untargeted/Untargeted_QCed_dat_long.csv"
matches.file <- "Intermediates/Untargeted/targeted_untargeted_MF_matches.csv"

###
dat <- read_csv(long.dat.file)

####get just 1 value for each triplicate
dat.1 <- dat %>%
  group_by(MF, method, Station) %>%
  mutate(mean_qc_area = mean(QCd_area)) %>%
  ungroup() %>%
  select(MF, Fraction, method, Retention_time, mz, Station, mean_qc_area) %>%
  mutate(Station = str_replace_all(Station, "A", "Aloha")) %>%
  unique() %>%
  mutate("method_station" = str_c(method, Station, sep = "_")) 

###bring in targeted matches to untargeted compounds
match.dat <- read_csv(matches.file)

dat.2 <- left_join(dat.1, match.dat)


dat.3 <- pivot_wider(dat.2, id_cols = c("MF", "Fraction", "mz", "Retention_time", "Compound.ID"), 
                     names_from = method_station, 
                     values_from = mean_qc_area) %>%
  arrange(mz) %>%
  select(MF, Fraction, mz, Retention_time, CatEx_Aloha, CatEx_PS, PPL_Aloha, PPL_PS, Compound.ID) %>%
  mutate(CatEx_Aloha = print(formatC(signif(CatEx_Aloha,digits=5), digits=5,format="d", flag="#"))) %>%
  mutate(CatEx_PS = print(formatC(signif(CatEx_PS,digits=5), digits=5,format="d", flag="#"))) %>%
  mutate(PPL_Aloha= print(formatC(signif(PPL_Aloha,digits=5), digits=5,format="d", flag="#"))) %>%
  mutate(PPL_PS= print(formatC(signif(PPL_PS,digits=5), digits=5,format="d", flag="#"))) %>%
  mutate(CatEx_Aloha = str_replace_all(CatEx_Aloha, "NA", "ND"),
         CatEx_PS = str_replace_all(CatEx_PS, "NA", "ND"),
         PPL_Aloha = str_replace_all(PPL_Aloha, "NA", "ND"),
         PPL_PS = str_replace_all(PPL_PS, "NA", "ND")) %>%
  mutate(Compound.ID = replace_na(Compound.ID, "unknown")) %>%
  rename("Mass Feature" = MF,
         "Analysis Fraction" = Fraction, 
         "m/z" = mz,
         "Retention Time (min)" = Retention_time,
         "CX-SPE - Aloha" = CatEx_Aloha,
         "CX-SPE - PS" = CatEx_PS,
         "PPL - Aloha" = PPL_Aloha,
         "PPL - PS" = PPL_PS,
         "Tentative Compound ID" = Compound.ID) 

#export
write_csv(dat.3, file = "Tables/Output/Untargeted_MFs_supptable.csv")
  
  
  
  
  
  
  
  
  