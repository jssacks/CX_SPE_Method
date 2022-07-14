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
library(ggpubr)


###Define inputs:
CX.file <- "Intermediates/Analytical_Validation/AV_CX_HQ_dat.csv"
PPL.file <- "Intermediates/Analytical_Validation/AV_PPL_HQ_dat.csv"
Std.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
ANOVA.file <- "Intermediates/Analytical_Validation/AV_CX_ANOVA_Results.csv"



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
  select(MF, Column, RT, m.z)  %>%
  mutate(MF = str_replace_all(.$MF, "Isoleucine", "(Iso)leucine")) %>%
  filter(!MF == "leucine")

All.dat.2 <- left_join(All.dat, Std.dat) %>%
  mutate(Fraction = case_when(
    .$Fraction == "Pos" ~ "HILIC Pos",
    .$Fraction == "Neg" ~ "HILIC Neg",
    TRUE ~ Fraction
  )) %>%
  mutate(method = case_when(
    .$method == "CatEx" ~ "CX-SPE",
    .$method == "PPL" ~ "PPL-SPE"
  )) %>%
  mutate("Fraction" = as.factor(Fraction))
All.dat.2$Fraction <- relevel(All.dat.2$Fraction, "HILIC Pos") 


#### MZ + RT + EE Plot

#figure data
Figure.dat <- All.dat.2 %>%
  mutate(MF = str_replace_all(.$MF, "Trimethylamine N-oxide", "TMAO"))

#figure palette
pal <- wes_palette("Zissou1", 100, type = "continuous")

#figure lables
label.comps <- tibble(MF = c("TMAO",
                               "Homarine",
                             "Serine",
                             "Alanine",
                             "Glutamic acid"
                           #  "Aspartic acid",
                           #  "Proline",
                           #  "DMSP",
                           #  "Betaine"
                             ))
Figure.labels <- inner_join(Figure.dat, label.comps)

#Make Figure
EE.fig <- ggplot(data = Figure.dat, aes(x = RT, y = m.z)) +
  geom_point(size = 2, alpha = 0.7, aes(color = Sample.Mean.EE, shape = Fraction)) +
  facet_grid(method~Column, scales = "free_x", space = "free_x") +
  scale_x_continuous(breaks = seq(0,20,3)) +
  theme_bw()+
  theme(panel.border = element_rect(colour = "black", size = 0.5), panel.grid = element_line(size = 0.1), axis.line = element_line(size = 0.5)) +
  theme(strip.background = element_rect(fill = "white", colour = "white")) +
  theme(legend.position = "bottom", legend.box.spacing = unit(0.2, "cm"), legend.title = element_text(size = 8), legend.text = element_text(size = 5.5)) +
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 6)) +
  theme(strip.text = element_text(size = 9)) + 
  guides(shape = "none") + 
  xlab("Retention Time (min)") +
  ylab("m/z") +
  labs(color = "Extraction Efficiency (%)") +
  geom_text_repel(data = Figure.labels, aes(x = RT, y = m.z, label = MF), 
                  box.padding = 1.3, min.segment.length = 0, ylim = 400, 
                  size = 2, 
                  segment.size = 0.2, 
                  segment.curvature = -1e-20) +
  #scale_color_gradient2(low = muted("lightskyblue1"), high = "gold2")
  scale_color_gradientn(colors = pal, trans = "sqrt", breaks = seq(0, 180, 20), guide = guide_colourbar(nbin = 100, barwidth = 12.5, barheight = 0.75)) #+
#  guides(color = guide_colourbar(nbin = 2, ticks = 8))
EE.fig
ggsave(filename = "Figures/Outputs/EE_Plot.pdf", width = 5, height = 4)



####ANOVA figures
ANOVA.dat <- read_csv(ANOVA.file)
parcor.dat <- ANOVA.dat %>%
  rowwise() %>%
  mutate(EE.norm = as.numeric(Sample.Mean.EE)/as.numeric(Overall.Mean.EE)) %>%
  mutate(RF.norm = as.numeric(Sample.Mean.RF)/as.numeric(Overall.Mean.RF)) %>%
  mutate(sample = as.factor(sample)) %>%
  select(MF, Fraction, sample, Sample.Mean.EE, EE.norm, Signif.EE, RF.norm, Signif.RF) %>%
  unique()

EE.parcor.plot <- ggplot(data = parcor.dat, aes(x = sample, y = EE.norm, group = MF, color = Signif.EE)) +
  geom_line(aes(alpha = Signif.EE)) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  scale_color_manual(values = c("black", "blue")) +
  theme_bw() +
  labs(color = "Signficant\nDifferences\nAcross\nSamples", alpha = "Signficant\nDifferences\nAcross\nSamples") +
  ylab("Normalized Extraction Efficiency") +
  geom_hline(aes(yintercept = 1), color = "red") +
  facet_wrap(.~Fraction, scales = "free_x")
EE.parcor.plot


RF.parcor.plot <- ggplot(data = parcor.dat, aes(x = sample, y = RF.norm, group = MF, color = Signif.RF)) +
  geom_line(aes(alpha = Signif.RF)) +
  scale_alpha_discrete(range = c(0.3, 1)) +
  scale_color_manual(values = c("black", "blue")) +
  theme_bw() +
  labs(color = "Signficant\nDifferences\nAcross\nSamples", alpha = "Signficant\nDifferences\nAcross\nSamples") +
  ylab("Normalized Response Factor") +
  geom_hline(aes(yintercept = 1), color = "red") +
  facet_wrap(.~Fraction, scales = "free_x")
RF.parcor.plot

#####
ggarrange(EE.parcor.plot, RF.parcor.plot, labels = c("A", "B"), common.legend = TRUE, nrow = 2, heights = 5, align = "v")



