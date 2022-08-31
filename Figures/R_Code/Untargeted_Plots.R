
#
#
#
#
#This script...
#
#
#
#
######
library(tidyverse)
library(readr)
library(ggpubr)
library(ggvenn)

###input
long.dat.file <- "Intermediates/Untargeted/Untargeted_QCed_dat_long.csv"

###load in data
dat <- read_csv(long.dat.file)

####get just 1 value for each triplicate
dat.1 <- dat %>%
  select(MF, Fraction, method, Fraction, Retention_time, mz, Station) %>%
  unique()

###Organize, regroup and tidy data
dat.2 <- dat.1 %>%
  mutate(chrom.col = case_when(str_detect(.$Fraction, "HILIC") ~ "HILIC",
                               str_detect(.$Fraction, "RP") ~ "RP")) %>%
  mutate(log.rt = log10(Retention_time),
         log.mz = log10(mz)) %>%
  mutate(method = str_replace(.$method, "CatEx", "CX-SPE")) %>%
  mutate(method = str_replace(.$method, "PPL", "PPL-SPE")) %>%
  mutate(Station = str_replace_all(.$Station, "A", "ALOHA"))

###Retention Time jitter and box plots of untargeted features
rt.plot <- ggplot(data = dat.2, aes(x = Retention_time, y = method)) +
  geom_jitter(alpha = 0.4, aes(), height = 0.20, size = 0.15) +
  geom_boxplot(alpha = 0.2, width = 0.5, size = 0.4, aes(color = method), outlier.size = 0.4) +
  scale_color_manual(values = c("darkred", "darkblue")) +
  facet_grid(Station~chrom.col, scales = "free") +
  theme_classic() + 
  xlab("Retention Time (min)") +
  theme(strip.background = element_rect(color = "white"),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9),
        axis.title = element_text(size = 7), 
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))
rt.plot


###m/z jitter and box plots of untargeted features
mz.plot <- ggplot(data = dat.2, aes(x = mz, y = method)) +
  geom_jitter(alpha = 0.4, aes(), height = 0.20, size = 0.15) +
  geom_boxplot(alpha = 0.2, width = 0.5, size = 0.4, aes(color = method), outlier.size = 0.4) +
  scale_color_manual(values = c("darkred", "darkblue")) +
  facet_grid(Station~chrom.col, scales = "free") +
  theme_classic() +
  xlab("m/z") +
  theme(axis.title.x = element_text(face = "italic"),
        strip.background = element_rect(color = "white"),
        strip.text.x = element_text(size = 9),
        strip.text.y = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 7), 
        axis.text = element_text(size = 6),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))
mz.plot

####Make final combined plot
ggarrange(rt.plot, mz.plot, labels = c("A", "B"), legend = "none", common.legend = TRUE, nrow = 2, heights = 5)
#save
ggsave(filename = "Figures/Outputs/untargeted_maintext_fig.pdf", height = 5, width = 5)

ggsave(filename = "Figures/Outputs/untargeted_maintext_fig.jpg", height = 5, width = 5, dpi = 300)


###Summary VENN DIAGRAM for Supplemental Figures

####Organize RP data for RP Venn Diagram 
dat.1.wide.RP <- dat.1 %>%
  unite(c("method", "Station"), col = method_station) %>%
  mutate(present = 1) %>%
  pivot_wider(id_cols = c(MF, Fraction), names_from = method_station, values_from = present) %>%
  replace(is.na(.), 0) %>%
  filter(str_detect(.$MF, "RP"))

###get MF IDs present in each of the 4 samples
dat.CXC.A <- dat.1.wide.RP %>%
  filter(CatEx_A == 1) %>%
  select(MF) %>%
  rename(`CX ALOHA` = MF)

dat.CXC.PS <- dat.1.wide.RP %>%
  filter(CatEx_PS == 1) %>%
  select(MF) %>%
  rename(`CX PS` = MF)

dat.PPL.A <- dat.1.wide.RP %>%
  filter(PPL_A == 1) %>%
  select(MF) %>%
  rename(`PPL ALOHA` = MF)

dat.PPL.PS <- dat.1.wide.RP %>%
  filter(PPL_PS == 1) %>%
  select(MF) %>%
  rename(`PPL PS` = MF)

#create list of list of IDs present in each method/sample pair
dat.RP <- as.list(c(dat.CXC.A, dat.CXC.PS, dat.PPL.A, dat.PPL.PS))

#Make Venn diagram
v.rp<- ggvenn(dat.RP,
              stroke_size = 0.25,
              set_name_size = 2,
              fill_alpha = 0.25,
              fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
              text_size = 1.75,
              show_percentage = TRUE)
v.rp




####Work on data set for HILIC Figure
dat.1.wide.HILIC <- dat.1 %>%
  unite(c("method", "Station"), col = method_station) %>%
  mutate(present = 1) %>%
  pivot_wider(id_cols = c(MF, Fraction), names_from = method_station, values_from = present) %>%
  replace(is.na(.), 0) %>%
  filter(!str_detect(.$MF, "RP"))

###get MF IDs present in each of the 4 samples
dat.CXC.A <- dat.1.wide.HILIC %>%
  filter(CatEx_A == 1) %>%
  select(MF) %>%
  rename(`CX ALOHA` = MF)

dat.CXC.PS <- dat.1.wide.HILIC %>%
  filter(CatEx_PS == 1) %>%
  select(MF) %>%
  rename(`CX PS` = MF)

dat.PPL.A <- dat.1.wide.HILIC %>%
  filter(PPL_A == 1) %>%
  select(MF) %>%
  rename(`PPL ALOHA` = MF)

dat.PPL.PS <- dat.1.wide.HILIC %>%
  filter(PPL_PS == 1) %>%
  select(MF) %>%
  rename(`PPL PS` = MF)

#create list of list of HILIC features in each sample/method pair
dat.HILIC <- as.list(c(dat.CXC.A, dat.CXC.PS, dat.PPL.A, dat.PPL.PS))

###plot hilic venn diagram 
v.hilic<- ggvenn(dat.HILIC,
                 stroke_size = 0.25,
                 set_name_size = 2.5,
                 fill_alpha = 0.25,
                 fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
                 text_size = 1.75,
                 show_percentage = TRUE) 


####Combine RP and HILIC Venn Diagrams into a single plot

ggarrange(v.hilic, v.rp, labels = c("A", "B"), nrow = 2)

ggsave(filename = "Figures/Outputs/Untargeted_Venn.pdf", height = 6, width = 5)








