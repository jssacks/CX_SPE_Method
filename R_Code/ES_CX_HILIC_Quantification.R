

#
#
#This script....


#####
library(readr)
library(tidyverse)


###Define inputs
norm.dat.file <- "Intermediates/Environmental_Samples/ES_CatEx_HILIC_BMISed_BlkQCed_dat.csv"
rf.file <- "Intermediates/Environmental_Samples/HILIC_RFs_RFratios.csv"
EE.file <- "Intermediates/Analytical_Validation/AV_CX_PPL_HQ_combined.csv"
Std.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
IS.names.file <- "Meta_Data/Ingalls_Standards/CXSPE_IS_List_JSS.csv"
HILIC.raw.dat <- "Intermediates/Environmental_Samples/ES_HILIC_targeted_combined_QCflags.csv"
Blk.LOD.dat <- "Intermediates/Environmental_Samples/ES_CX_HILIC_Blk_LOD_signal_values.csv"
IS.raw.file.pos <- "Intermediates/Environmental_Samples/ES_HILIC_pos_IS.csv"
IS.raw.file.neg <- "Intermediates/Environmental_Samples/ES_HILIC_neg_IS.csv"

###HILIC Dat
hilic.dat <- read_csv(norm.dat.file) %>%
  mutate(MF = str_replace_all(.$MF, "butyric_acid_Neg", "butyric acid_Neg")) %>%
  separate(MF, into = c("Compound", "ion_mode"), sep = "_") %>%
  mutate(z = ifelse(ion_mode == "Neg", -1, 1)) %>%
  mutate(sample = str_replace_all(.$samp, "Aloha", "A")) %>%
  mutate(Compound = str_replace_all(.$Compound, "Isoleucine", "(Iso)leucine")) %>%
  filter(!Compound == "leucine")

###Join HILIC Dat and RFs and RF ratios
hilic.rfs <- read_csv(rf.file) %>%
  mutate(Compound = str_replace_all(.$Compound, "Isoleucine", "(Iso)leucine")) %>%
  filter(!Compound == "leucine")


dat <- left_join(hilic.dat, hilic.rfs, by = c("Compound")) %>%
  filter(!is.na(RF))

###Calculate environmental concentration using RFs
dat.conc <- dat %>%
  rowwise() %>%
  mutate(RF = as.numeric(RF)) %>%
  mutate(nmol.conc = Adjusted_Area/RF/RFratio*10^-6*400/(40*10^-3)*1000)

###Adjust values using extraction efficiencies and filter by unretained compounds
EE.hilic <- read_csv(EE.file) %>%
  select(MF, method, Fraction, sample, Overall.Mean.EE, Overall.SD.EE, 
         Sample.Mean.EE, Sample.SD.EE) %>%
  filter(method == "CX-SPE") %>%
  filter(!Fraction == "RP") %>%
  rename("Compound" = MF) %>%
  mutate(Compound = str_replace_all(.$Compound, "Isoleucine", "(Iso)leucine")) %>%
  filter(!Compound == "leucine")

dat.conc.EE <- left_join(dat.conc, EE.hilic, by = c("Compound")) %>%
  filter(!is.na(Overall.Mean.EE)) %>%
  mutate(EE.adjust.conc = nmol.conc/(Overall.Mean.EE/100)) %>%
  unique()

####Calculate LODs in Concentration Space 
blk.lod.dat <- read_csv(Blk.LOD.dat) %>%
  mutate(MF = str_replace_all(.$MF, "butyric_acid_Neg", "butyric acid_Neg")) %>%
  separate(MF, into = c("Compound", "ion_mode"), sep = "_") %>%
  mutate(z = ifelse(ion_mode == "Neg", -1, 1)) %>%
  mutate(LOQ.test = Blk.Av + (10 * (Blk.sd/sqrt(15)))) %>%
  mutate(Compound = str_replace_all(.$Compound, "Isoleucine", "(Iso)leucine")) %>%
  filter(!Compound == "leucine")

lod.dat.2 <- left_join(blk.lod.dat, hilic.rfs, by = c("Compound", "z")) %>%
  filter(!is.na(RF))
  

lod.conc <- lod.dat.2 %>%
  rowwise() %>%
  mutate(RF = as.numeric(RF)) %>%
  mutate(lod.nmol.conc = Blk.LD/RF/RFratio*10^-6*400/(40*10^-3)*1000,
         loq.nmol.conc = LOQ.test/RF/RFratio*10^-6*400/(40*10^-3)*1000)

lod.conc.EE <-  left_join(lod.conc, EE.hilic, by = c("Compound")) %>%
  filter(!is.na(Overall.Mean.EE)) %>%
  mutate(EE.adjust.lod = lod.nmol.conc/(Overall.Mean.EE/100),
         EE.adjust.loq = loq.nmol.conc/(Overall.Mean.EE/100)) %>%
  select(Compound, z, EE.adjust.lod, EE.adjust.loq) %>%
  unique()

####Calculate better concentrations using Internal standards
IS_names <- read_csv(IS.names.file) %>%
  rename(Compound = Match.Old,
         IS = IS.Name.Old)

#bring in raw IS data and combine with IS matching
IS.dat.pos <- read_csv(IS.raw.file.pos)
IS.dat.neg <- read_csv(IS.raw.file.neg)
IS.dat.full <- rbind(IS.dat.pos, IS.dat.neg) %>%
  rename(IS = Compound,
         IS_Area = Area)
IS.dat.named <- left_join(IS.dat.full, IS_names)

###Spike Before:
IS.Spike.Before <- IS.dat.named %>%
  filter(Spike.Fraction == "Spike_Before")

raw.dat <- read_csv(HILIC.raw.dat) %>%
  separate(Compound, into = c("Compound", "ion_mode")) %>%
  mutate(z = ifelse(ion_mode == "Neg", -1, 1)) %>%
  select(Rep, Compound, z, Area)

SBe.Matched <- left_join(IS.Spike.Before, raw.dat, by = c("Compound", "Rep")) %>%
  filter(!str_detect(.$Rep, "Std")) %>%
  filter(!str_detect(.$Rep, "Blk")) %>%
  filter(!str_detect(.$Rep, "_C_")) %>%
  filter(!str_detect(.$Rep, "PPL")) %>%
  filter(!str_detect(.$Rep, "Poo")) %>%
  select(Rep, Compound, Area, IS_Area, Conc.in.vial.uM) %>%
  mutate(Nmol.in.vial_IS = Area/IS_Area*Conc.in.vial.uM*1000,
         Nmol.in.Samp_IS = Nmol.in.vial_IS*10^-6*400/(40*10^-3))

###Manually verified quality of IS calcs, removing bad IS calcs
SBe.Checked <- SBe.Matched %>%
 # filter(!Compound == "Homarine") %>%
  filter(!Compound == "Histidine") %>%
  filter(!Compound == "Methionine") %>%
  filter(!Compound == "Sucrose") %>%
  filter(!Compound == "Valine") %>%
  filter(!Compound == "Guanine")

SBe.add <- SBe.Checked %>%
  rename(EE.adjust.conc = Nmol.in.Samp_IS) %>%
  select(Compound, Rep, EE.adjust.conc)

####Spike After Internal Standards 
IS.Spike.After <- IS.dat.named %>%
  filter(Spike.Fraction == "Spike_After")

SAf.Matched <- left_join(IS.Spike.After, raw.dat, by = c("Compound", "Rep")) %>%
  filter(!str_detect(.$Rep, "Std")) %>%
  filter(!str_detect(.$Rep, "Blk")) %>%
  filter(!str_detect(.$Rep, "_C_")) %>%
  filter(!str_detect(.$Rep, "PPL")) %>%
  filter(!str_detect(.$Rep, "Poo")) %>%
  select(Rep, Compound, Area, IS_Area, Conc.in.vial.uM) %>%
  mutate(Nmol.in.vial_IS = Area/IS_Area*Conc.in.vial.uM*1000,
         Nmol.in.Samp_IS = Nmol.in.vial_IS*10^-6*400/(40*10^-3))

SAf.EE <- left_join(SAf.Matched, EE.hilic) %>%
  filter(!is.na(Overall.Mean.EE)) %>%
  mutate(EE.adjust.conc.IS = Nmol.in.Samp_IS/(Overall.Mean.EE/100)) %>%
  select(Rep, Compound, EE.adjust.conc.IS) %>%
  unique()

SAf.EE.add <- SAf.EE %>%
  rename(EE.adjust.conc = EE.adjust.conc.IS)

###Combine IS normalized Data, compare to LODs, and prepare it for adding it back into the larger data set

#combine SAf and SBe IS
IS.adjus.dat <- rbind(SAf.EE.add, SBe.add) %>%
  separate(Rep, 
           c("runDate",
             "type","samp","replicate"),"_", remove = FALSE)

###remove calculated Concentrations below LOD
LOD.IS.adjus.dat <- left_join(IS.adjus.dat, lod.conc.EE) %>%
  mutate(lod.flag = ifelse(EE.adjust.conc <= EE.adjust.lod, 1, 0),
         loq.flag = ifelse(EE.adjust.conc <= EE.adjust.loq, 1, 0))

LOD.IS.QC <- LOD.IS.adjus.dat %>%
  group_by(Compound, samp) %>%
  summarise(lod.flag.sum = sum(lod.flag)) %>%
  filter(!lod.flag.sum > 0) %>%
  mutate(lod.IS.QC = "ok") %>%
  select(Compound, samp, lod.IS.QC)

IS.dat.QCed <- left_join(LOD.IS.adjus.dat, LOD.IS.QC) %>%
  filter(lod.IS.QC == "ok") %>%
  select(Compound, Rep, samp, EE.adjust.conc) %>%
  rename(sample = samp) %>%
  mutate(sample = str_replace_all(.$sample, "Aloha", "A")) 


####Add in IS normalized dat
no.IS.dat <- dat.conc.EE %>%
  select(Compound, Rep, samp, EE.adjust.conc) %>%
  rename(sample = samp)

###
Remove.dat <- no.IS.dat %>%
  filter(!Compound %in% IS.dat.QCed$Compound)

##
final.conc.dat <- rbind(Remove.dat, IS.dat.QCed)

#####Get mols C and N per mol Compound from empirical formula in standard sheet
std.info <- read_csv(Std.file)

C.N.std.info <- std.info %>%
  filter(Priority = TRUE) %>%
  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
                    str_extract(Emperical.Formula, "^C\\d"), 
                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
                    1, str_extract(Emperical.Formula, "N\\d")))%>%
  mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
  select(Compound.Name_old, C, N) %>%
  rename(Compound = Compound.Name_old) %>%
  mutate(Compound = str_replace_all(.$Compound, "Isoleucine", "(Iso)leucine")) %>%
  filter(!Compound == "leucine")

###Conc Data with C and N mol space data
final.dat <- left_join(final.conc.dat, C.N.std.info) %>%
  rowwise() %>%
  mutate(Nmol.C = C*EE.adjust.conc,
         Nmol.N = N*EE.adjust.conc) %>%
  filter(!sample == "TruePoo") %>%
  unique()

###Write Environmental Concentraions and LOD Concentrations to a csv:
#Enviro.Concs:
write_csv(final.dat, file = "Intermediates/Environmental_Samples/ES_CX_HILIC_Concentrations.csv")

#LODs:
write_csv(lod.conc.EE, file = "Intermediates/Environmental_Samples/ES_CX_Blk_LOD_Concentrations.csv")

