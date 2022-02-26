#
#
#
#
#
#
#packages
library(tidyverse)
library(readr)
library(wesanderson)
library(ggrepel)


###Define inputs:
CX.file <- "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv"
PPL.file <- "Intermediates/Analytical_Validation/AV_PPL_HQ_dat.csv"
Std.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"


#####Bring in and join CX and PPL dat
CX.dat <- read_csv(CX.file) %>%
  mutate("method" = "CatEx") %>%
  select(MF, Fraction, Sample.Mean.EE, method)
PPL.dat <- read_csv(PPL.file) %>%
  mutate("method" = "PPL")  %>%
  select(MF, Fraction, Sample.Mean.EE, method)

All.dat <- rbind(CX.dat, PPL.dat)
  
###Add in compound specific information from standards sheet (mz and RT)
Std.dat <- read_csv(Std.file)  %>%
  mutate(MF = Compound.Name_old,
         Column = Column,
         RT = RT..min.) %>%
  mutate(m.z = as.numeric(m.z)) %>%
  filter(Priority == TRUE) %>%
  select(MF, Column, RT, m.z) 

All.dat.2 <- left_join(All.dat, Std.dat) %>%
  mutate("Fraction" = as.factor(Fraction))
All.dat.2$Fraction <- relevel(All.dat.2$Fraction, "Pos")


#### MZ + RT + EE Plot

#figure data
Figure.dat <- All.dat.2

#figure palette
pal <- wes_palette("Zissou1", 100, type = "continuous")

#figure lables
label.comps <- tibble(MF = c("Gonyol",
                               "Homarine",
                             "Ectoine",
                             #   "Hydroxyectoine",
                             "Phenylalanine",
                             "Domoic Acid",
                             "Vitamin B2",
                             "Thymine"))
Figure.labels <- inner_join(Figure.dat, label.comps)

#Make Figure
EE.fig <- ggplot(data = Figure.dat, aes(x = RT, y = m.z)) +
  geom_point(size = 4.0, alpha = 0.8, aes(color = Sample.Mean.EE, shape = Fraction)) +
  facet_grid(method~Column, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0,20,3)) +
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", size = 1)) +
  theme(strip.background = element_rect(fill = "white", colour = "white")) +
  xlab("Retention Time (min)") +
  ylab("m/z") +
  labs(color = "Extraction Efficiency (%)") +
  geom_text_repel(data = Figure.labels, aes(x = RT, y = m.z, label = MF), 
                  box.padding = 1.5, min.segment.length = 0, ylim = 500,
                  arrow = arrow(length = unit(0.015, "npc")),
                  segment.curvature = -1e-20) +
  #scale_color_gradient2(low = muted("lightskyblue1"), high = "gold2")
  scale_color_gradientn(colors = pal, trans = "sqrt", breaks = seq(0, 180, 20), guide = guide_colourbar(nbin = 100, barwidth = 1.5, barheight = 15)) #+
#  guides(color = guide_colourbar(nbin = 2, ticks = 8))
EE.fig










