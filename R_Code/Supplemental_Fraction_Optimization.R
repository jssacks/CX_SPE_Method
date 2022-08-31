#load packages
library(tidyverse)
library(readr)
library(stringr)
library(viridis)
#load functions
source("R_Code/CXSPE_Functions.R")


###Define data files
pos.file <- "Raw_Data/Optimization_Test/Cat_SPE_spike_results_POS_200113.csv"
neg.file <- "Raw_Data/Optimization_Test/Cat_SPE_spike_results_Neg_200115.csv"




####load in data
pos.dat <- csv_parse(pos.file) %>%
  mutate(mode = "Pos")
neg.dat <- csv_parse(neg.file) %>%
  mutate(mode = "Neg")


##combine data, get rid of negative mode because CX-SPE does not really work for these compounds,
# and kick out compounds that are not observed.
full.dat <- rbind(pos.dat, neg.dat) %>%
  mutate(Rep = str_remove(.$Rep, "200106_")) %>%
  filter(!str_detect(.$Rep, "Blk")) %>%
  mutate(Fraction = as.numeric(str_remove(.$Rep, "Elut_"))) %>%
  filter(!str_detect(.$Compound, ",")) %>%
  filter(!mode == "Neg") %>% 
  mutate(Area = replace_na(.$Area, 0)) %>%
  group_by(Compound) %>%
  mutate(mean.area = mean(Area),
         area.total = sum(Area),
         area.perc = (Area/area.total)*100) %>%
  filter(!mean.area <= 40000)

ggplot(full.dat, aes(x = factor(Fraction), y = Compound, fill = area.perc))  +
  geom_tile(colour="black", size=0.1) +
  scale_fill_viridis("Percentage of \nTotal Peak Area (%)") +
  xlab("Fraction") +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 6))
  

ggplot(full.dat, aes(x = Fraction, y = area.perc)) +
  geom_point(alpha = 0.4) +
  geom_boxplot(aes(group = Rep), alpha = 0.1)



sum.dat <- full.dat %>%
  ungroup() %>%
  group_by(Fraction) %>%
  summarise(Average.Perc.Area = mean(area.perc))



















