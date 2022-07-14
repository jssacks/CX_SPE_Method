###Estimate the number of targeted compounds found in the untargeted dataset 


library(tidyverse)
library(readr)
library(ggpubr)


###Define thresholds

#RT threshold (min)
RT.error <- 0.05

#PPM threshold (ppm)
PPM.error <- 5




###input
#untargeted data file
long.dat.file <- "Intermediates/Untargeted/Untargeted_QCed_dat_long.csv"

#targeted compound details
RP.tlist.file <- "Meta_Data/Transition_Lists/RP_TransitionList.csv"
HPos.tlist.file<- "Meta_Data/Transition_Lists/HILIC_Pos_TransitionList.csv"
HNeg.tlist.file <- "Meta_Data/Transition_Lists/HILIC_Neg_TransitionList.csv"

#targeted data files:


#standard file
Std.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"


#load functions
source("R_Code/CXSPE_Functions.R")


#Define inputs
#HILIC positive skyline ES files for CX-SPE
sky.file.pos <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_Pos_G3.csv"

#HILIC negative skyline ES files for CX-SPE
sky.file.neg <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_Neg_G3.csv"

#RP skyline ES file for CX-SPE
sky.file.rp <- "Raw_Data/Skyline_Output/Environmental_Samples/CXSPE_ES_RP_G1.csv"

#sample list
smp.list <- "Meta_Data/Sample_Lists/CXSPE_ES_HILIC_SampList.csv"


#################### ___________________________________________
#####define functions to load and parse files
sky_read <- function(file) {
  output <- read_csv(file) %>%
    mutate(Rep =`Replicate Name`,
           Compound = `Precursor Ion Name`,
           RT = as.numeric(`Retention Time`),
           Area = as.numeric(`Area`),
           Background = as.numeric(`Background`),
           Height = as.numeric(`Height`),
           Mass_Error_PPM = as.numeric(`Mass Error PPM`)
    ) %>%
    select(Rep, Compound, RT, Area, Background, Height, Mass_Error_PPM)
  print(output)
}



tl_read <- function(file) {
  output <- read_csv(file, col_names = FALSE) %>%
    mutate(Mass = as.numeric(X2),
           Compound = X4) %>%
    select(Compound, Mass)
}

###get new and old names from standard sheet
Std.name.key <- read_csv(Std.file)  %>%
  select(Compound.Name, Compound.Name_old) %>% 
  rename("Compound" = Compound.Name_old) %>%
  unique()


########load in files
HPos.dat <- sky_read(sky.file.pos)
Hpos.dat.fixednames <- left_join(HPos.dat, Std.name.key) %>%
  select(-Compound) %>%
  rename("Compound" = Compound.Name)
HPos.tl <- tl_read(HPos.tlist.file) %>%
  mutate(Compound = str_replace_all(.$Compound, "_", ","))

HNeg.dat <- sky_read(sky.file.neg)
HNeg.dat.fixednames <- left_join(HNeg.dat, Std.name.key) %>%
  select(-Compound) %>%
  rename("Compound" = Compound.Name)
HNeg.tl <- tl_read(HNeg.tlist.file) %>%
  mutate(Compound = str_replace_all(.$Compound, "_", ","))

RP.dat <- sky_read(sky.file.rp)
RP.dat.fixednames <- left_join(RP.dat, Std.name.key) %>%
  select(-Compound) %>%
  rename("Compound" = Compound.Name)
RP.tl <- tl_read(RP.tlist.file) %>%
  mutate(Compound = str_replace_all(.$Compound, "_", ","))

#####combine pooled data and transition list information
## calculate batch specific exact mass by incorporating batch specific mass defect

####HPos
HPos.dat.em <- left_join(Hpos.dat.fixednames, HPos.tl) %>%
  filter(str_detect(.$Rep, "Poo")) %>%
  mutate(Exact_Mass = Mass + (Mass_Error_PPM/10^6)*Mass) %>%
  group_by(Compound) %>%
  summarize(Mean.Exact.Mass = mean(Exact_Mass), 
            Mean.RT = mean(RT)) %>%
  rename("Exact_Mass" = Mean.Exact.Mass,
         "RT" = Mean.RT) %>%
  mutate(Fraction = "HILICPos")  %>%
  filter(!Exact_Mass == "NaN") %>%
  filter(!RT == "NaN")

######HNeg
HNeg.dat.em <- left_join(HNeg.dat.fixednames, HNeg.tl) %>%
  filter(str_detect(.$Rep, "Poo")) %>%
  mutate(Exact_Mass = Mass + (Mass_Error_PPM/10^6)*Mass) %>%
  group_by(Compound) %>%
  summarize(Mean.Exact.Mass = mean(Exact_Mass, na.rm = TRUE), 
            Mean.RT = mean(RT, na.rm = TRUE)) %>%
  rename("Exact_Mass" = Mean.Exact.Mass,
         "RT" = Mean.RT) %>%
  mutate(Fraction = "HILICNeg")  %>%
  filter(!Exact_Mass == "NaN") %>%
  filter(!RT == "NaN")

####RP
RP.dat.em <- left_join(RP.dat.fixednames, RP.tl) %>%
  filter(str_detect(.$Rep, "Poo")) %>%
  mutate(Exact_Mass = Mass + (Mass_Error_PPM/10^6)*Mass) %>%
  group_by(Compound) %>%
  summarize(Mean.Exact.Mass = mean(Exact_Mass, na.rm = TRUE), 
            Mean.RT = mean(RT, na.rm = TRUE)) %>%
  rename("Exact_Mass" = Mean.Exact.Mass,
         "RT" = Mean.RT) %>%
  mutate(Fraction = "RP") %>%
  filter(!Exact_Mass == "NaN") %>%
  filter(!RT == "NaN")


##Create targeted dataset
dat.target <- rbind(HPos.dat.em, HNeg.dat.em, RP.dat.em) %>%
  mutate(Exact_Mass_low = Exact_Mass - (PPM.error/10^6)*Exact_Mass,
         Exact_Mass_high = Exact_Mass + (PPM.error/10^6)*Exact_Mass,
         RT_low = RT - RT.error,
         RT_high = RT + RT.error)

###Bring in untargeted data
untarg.dat <- read_csv(long.dat.file) 



####Find targeted compounds in untargeted data!

###Find Targeted Compound Function
find_targeted <- function(target.comp, Skyline.df, MSDIAL.df) {
  # 
  target.df <- Skyline.df %>%
    filter(Compound == target.comp)
  #
  untargeted.df <- MSDIAL.df %>%
    filter(Fraction == target.df$Fraction) %>%
    filter(between(mz, target.df$Exact_Mass_low, target.df$Exact_Mass_high)) %>%
    filter(between(Retention_time, target.df$RT_low, target.df$RT_high)) %>%
    mutate(Compound.ID = target.df$Compound)
  #
  print(untargeted.df)
}


#create ouput dataframe
found.targets <- data.frame(
  "Compound.ID" = character(),
  "Fraction" = character(),
  "MF" = character(),
  "SampID" = character(),
  "Retention_time" = numeric(),
  "mz" = numeric(),
  "QCd_area" = numeric(),
  "method" = character(),
  "Station" = character()
)

###loop through skyline masses to find data
sky.data.loop <- dat.target
msdial.data.loop <- untarg.dat

for (i in seq_along(sky.data.loop$Compound)) {
  comp <- sky.data.loop$Compound[i]
  msdial.out <- find_targeted(comp, sky.data.loop, msdial.data.loop)
  msdial.out.2 <- msdial.out %>%
    mutate("Compound" = comp) %>%
    select(Compound.ID, Fraction, MF, SampID, Retention_time, mz, QCd_area, method, Station) 
  found.targets <- full_join(found.targets, msdial.out.2)
}


###Look at results!
out.loop <- found.targets
target.comps <- out.loop %>%
  select(Compound.ID, MF, Fraction) %>%
  unique()

found.comps.summary <- target.comps %>%
  group_by(Fraction) %>%
  summarize(count = n())

write_csv(target.comps, file = "Intermediates/Untargeted/targeted_untargeted_MF_matches.csv")






