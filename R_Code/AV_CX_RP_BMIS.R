library(ggplot2)
library(tidyverse)


#Define all your inputs here
samp.key.file <- "Intermediates/Analytical_Validation/AV_CX_RP_Samp_List.csv"
is.dat.file <- "Intermediates/Analytical_Validation/AV_RP_IS.csv"
xcms.dat.file <- "Intermediates/Analytical_Validation/AV_RP_combined_raw.csv"
cut.off <- 0.2
cut.off2 <- 0.1


#Import sample key----
samp.key <- read_csv(samp.key.file) 

###Added in filter to remove PPL Samples
xcms.dat.2 <- read_csv(xcms.dat.file) %>%
  rename("MF" = Compound) %>%
  pivot_longer(-MF, names_to = "SampID", values_to = "Area") %>%
  filter(!str_detect(.$SampID, "PPL")) 

xcms.dat <- xcms.dat.2 %>%
  mutate(Area = as.numeric(Area))

#Import and clean up the Internal standard data----
is.dat <- read_csv(is.dat.file) %>%
  rename(SampID = `Rep`,
         MF = `Compound`) %>% select(SampID, MF, Area) %>%
  filter(!str_detect(.$SampID, "PPL"))


is.dat.full <- is.dat

#Read in sample key, add in injec_volume data from sample key----
samp.key.2 <- samp.key %>%
  filter(Rep %in% is.dat.full$SampID) %>%
  select(Rep, InjVol_norm) %>%
  filter(!is.na(InjVol_norm))%>%
  mutate(MF = "Inj_vol",
         Area = InjVol_norm,
         SampID = Rep) %>%
  select(SampID, MF, Area)

is.dat.full.with.samp <- rbind(is.dat.full, samp.key.2)

#Look at extraction replication of the Internal Standards----
IS_inspectPlot <- ggplot(is.dat.full.with.samp, aes(x=SampID, y=Area)) + 
  geom_bar(stat="identity") + 
  facet_wrap( ~MF, scales="free_y")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 5), 
        axis.text.y = element_text(size = 10),
        legend.position = "top",
        strip.text = element_text(size = 10))+
  ggtitle("IS Raw Areas")
IS_inspectPlot


#Edit data so names match, separate out the date so we can get individual IS.means for each run batch, remove any ISs that we don't think are good to use
is.dat.full.with.samp.edited <- is.dat.full.with.samp %>%
  filter(str_detect(MF, "L-Phenylalanine") |
           str_detect(MF, "Riboflavin") |
           str_detect(MF, "Pyridoxal") |
           str_detect(MF, "Inj_vol"))

###Change name of data frame
xcms.long <- xcms.dat



#Calculate mean values for each IS----
is.means <- is.dat.full.with.samp.edited %>% 
  left_join(samp.key %>%
              mutate(SampID = Rep) %>% 
              select(SampID), by = "SampID") %>%
  group_by(MF) %>%
  summarise(ave = mean(as.numeric(Area))) %>%
  mutate(ave = ifelse(MF == "Inj_vol", 1, ave))
#REMOVED SAMPLE.TYPE


#Normalize to each internal Standard----
binded <- rbind(is.dat.full.with.samp.edited, xcms.long) %>%
  left_join(samp.key %>%
              mutate(SampID = Rep) %>% 
              select(SampID), by = "SampID")
split.dat <- list()
for (i in 1:length(unique(is.dat.full.with.samp.edited$MF))){
  split.dat[[i]] <- binded %>% mutate(MIS = unique(is.dat.full.with.samp.edited$MF)[i]) %>%
    left_join(is.dat.full.with.samp.edited %>% 
                rename(MIS = MF, IS_Area = Area) %>% 
                select(MIS, SampID, IS_Area)) %>%
    left_join(is.means %>% 
                rename(MIS = MF)) %>%
    mutate(Adjusted_Area = Area/IS_Area*ave)
}
area.norm <- do.call(rbind, split.dat) %>% select(-IS_Area, -ave)


#Break Up the Names (Name structure must be:  Date_type_ID_replicate_anythingextraOK)----
area.norm.2 <- area.norm %>% separate(SampID, 
                                      c("runDate",
                                        "type","samp","replicate"),"_", remove = FALSE)

#Find the B-MIS for each MassFeature----
#Look only the Pooled samples, to get a lowest RSD of the pooled possible (RSD_ofPoo), 
#then choose which IS reduces the RSD the most (Poo.Picked.IS) 
poodat <- area.norm.2 %>%
  filter(type == "Poo")%>% 
  group_by(samp, MF, MIS) %>%
  summarise(RSD_ofPoo_IND = sd(Adjusted_Area, 
                               na.rm = TRUE)/mean(Adjusted_Area, na.rm = TRUE))%>%
  mutate(RSD_ofPoo_IND = ifelse(RSD_ofPoo_IND == "NaN", NA, RSD_ofPoo_IND)) %>%
  group_by(MF, MIS) %>%
  summarise(RSD_ofPoo =  mean(RSD_ofPoo_IND, na.rm = TRUE))

poodat <- poodat %>% left_join(poodat %>%
                                 group_by(MF) %>%
                                 summarise(Poo.Picked.IS = 
                                             unique(MIS)[which.min(RSD_ofPoo)][1]))

#Get the starting point of the RSD (Orig_RSD), calculate the change in the RSD, say if the MIS is acceptable----
poodat.2 <- left_join(poodat, poodat %>%
                        filter(MIS == "Inj_vol" ) %>%
                        mutate(Orig_RSD = RSD_ofPoo) %>%
                        select(-RSD_ofPoo, -MIS)) %>%
  mutate(del_RSD = (Orig_RSD - RSD_ofPoo)) %>%
  mutate(percentChange = del_RSD/Orig_RSD) %>%
  mutate(accept_MIS = (percentChange > cut.off & Orig_RSD > cut.off2)) 

#Change the BMIS to "Inj_vol" if the BMIS is not an acceptable -----
#Adds a column that has the BMIS, not just Poo.picked.IS
#Changes the finalBMIS to inject_volume if its no good
fixedpoodat <- poodat.2 %>%
  filter(MIS == Poo.Picked.IS)%>%
  mutate(FinalBMIS = ifelse((accept_MIS == "FALSE"), "Inj_vol", Poo.Picked.IS), 
         FinalRSD = RSD_ofPoo) 
newpoodat <- poodat.2 %>% left_join(fixedpoodat %>% select(MF, FinalBMIS)) %>%
  filter(MIS == FinalBMIS) %>%
  mutate(FinalRSD = RSD_ofPoo)
report.text <- newpoodat %>% filter(FinalBMIS != "Inj_vol")
QuickReport <- paste("% of MFs that picked a BMIS", 
                     length(report.text$MF) / length(newpoodat$MF), 
                     "RSD improvement cutoff", cut.off,
                     "RSD minimum cutoff", cut.off2,
                     sep = " ")
QuickReport

#Evaluate the results of your BMIS cutoff-----
IS_toISdat <- area.norm.2 %>%
  filter(MF %in% is.dat.full.with.samp.edited$MF) %>%
  select(MF, MIS, Adjusted_Area, type) %>%
  filter(type == "Smp") %>%
  group_by(MF, MIS) %>%
  summarise(RSD_ofSmp = sd(Adjusted_Area)/mean(Adjusted_Area)) %>%
  left_join(poodat %>% select(MF, MIS, RSD_ofPoo))

injectONlY_toPlot <- IS_toISdat %>%
  filter(MIS == "Inj_vol" ) 


ISTest_plot <- ggplot()+
  geom_point(dat = IS_toISdat, shape = 21, color = "black", size = 2,aes(x = RSD_ofPoo, y = RSD_ofSmp))+ 
  scale_fill_manual(values=c("white","dark gray"))+
  geom_point(dat = injectONlY_toPlot, aes(x = RSD_ofPoo, y = RSD_ofSmp), size = 3) +
  facet_wrap(~ MF)
ISTest_plot

#Get all the data back - and keep only the MF-MIS match set for the BMIS----
#Add a column to the longdat that has important information from the FullDat_fixed, 
#then only return data that is normalized via B-MIS normalization
BMIS_normalizedData <- newpoodat %>% select(MF, FinalBMIS, Orig_RSD, FinalRSD) %>%
  left_join(area.norm.2 %>% rename(FinalBMIS = MIS)) %>% unique() %>%
  filter(!MF %in% is.dat.full.with.samp.edited$MF)


BMIS_normalizedData.2 <- BMIS_normalizedData 
write_csv(BMIS_normalizedData.2, "Intermediates/Analytical_Validation/AV_CX_RP_BMISed_dat.csv")

BMISlist <- list(IS_inspectPlot, QuickReport, ISTest_plot, BMIS_normalizedData.2)
ISTest_plot
#Removes all intermediate variables :)
rm(list=setdiff(ls(), c("BMISlist")))
