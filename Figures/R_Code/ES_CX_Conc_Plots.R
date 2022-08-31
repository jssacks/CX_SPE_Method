library(tidyverse)
library(readr)
library(ggthemes)
library(palettetown)
library(ggsci)
library(jcolors)
library(LaCroixColoR)
library(ggpubr)

#install.packages("devtools")
#devtools::install_github("johannesbjork/LaCroixColoR")
####Define inputs
HILIC.conc.file <- "Intermediates/Environmental_Samples/ES_CX_HILIC_Concentrations.csv"
RP.conc.file <- "Intermediates/Environmental_Samples/ES_CXC_RP_Concentrations.csv"
#HILIC.LOD.file <- "Intermediates/Environmental_Samples/ES_Blk_LOD_Concentrations.csv"
#RP.LOD.file <- "Intermediates/Environmental_Samples/ES_CXC_RP_BLk_LOD_Concentrations.csv"
###
HILIC.conc <- read_csv(HILIC.conc.file) %>%
  mutate(sample = case_when(.$sample == "A" ~ "ALOHA",
                            .$sample == "Aloha" ~ "ALOHA",
                            TRUE ~ sample))

RP.conc <- read_csv(RP.conc.file) %>%
  mutate(sample = str_replace_all(.$samp, "A", "ALOHA")) %>%
  mutate(sample = str_replace_all(.$sample, "ALOHAloha", "ALOHA")) %>%
  select(-samp)

##
dat <- rbind(HILIC.conc, RP.conc) %>%
  filter(!Compound == "Salicylic Acid") %>%
  group_by(Compound, sample) %>%
  mutate(mean.conc = mean(EE.adjust.conc),
         sd.mean.conc = sd(EE.adjust.conc),
         mean.Nmol.C = mean(Nmol.C),
         sd.Nmol.C = sd(Nmol.C),
         mean.Nmol.N = mean(Nmol.N),
         sd.Nmol.N = sd(Nmol.N)) %>%
  select(-EE.adjust.conc, -Nmol.C, -Nmol.N, -Rep) %>%
  unique()

####Supplementary Table of Enviro Conc. 
dat.table <- dat %>%
  select(-C, -N)
write_csv(dat.table, path = "Intermediates/CX_Enviro_Conc_Supp_Table.csv")

dat.sum <- dat %>%  
  group_by(sample) %>%
  summarise(mean.conc = sum(mean.conc),
            Nmol.C = sum(mean.Nmol.C),
            Nmol.N = sum(mean.Nmol.N, na.rm = TRUE),
            count = n()) 
dat.sum
###get 12 most abundant compounds
ma.dat <- dat %>%
  group_by(Compound) %>%
  filter(!Compound == "Salicylic Acid") %>%
  filter(!Compound == "Ornithine") %>%
  summarize(max.abu = max(mean.Nmol.C)) %>%
  slice_max(max.abu, n = 11) %>%
  mutate(order = row_number()) %>%
  select(Compound, order)

####Summarize the mean.conc, nmol.C, and nmol.N of the less abundant compounds
dat.other <- full_join(dat, ma.dat) %>%
  filter(is.na(order)) %>%
  group_by(sample) %>%
  summarise(mean.conc = sum(mean.conc),
            mean.Nmol.C = sum(mean.Nmol.C),
            mean.Nmol.N = sum(mean.Nmol.N, na.rm = TRUE)) %>%
  mutate(Compound = "Other",
         order = 12)



###Put most abundant and other data together for figure  
ma.dat.fig <- full_join(dat, ma.dat) %>%
  filter(!is.na(order)) %>%
  select(Compound, sample, mean.conc, mean.Nmol.C, mean.Nmol.N, order)

dat.fig <- rbind(ma.dat.fig, dat.other) 

#Change Trimethylaammonium Propionate (TMAP) to just TMAP
dat.fig <- dat.fig %>%
  ungroup() %>%
  mutate(across(`Compound`, str_replace, "Trimethylammonium Propionate", "TMAP")) %>%
  mutate(across(`Compound`, str_remove, "\\(TMAP\\)")) %>%
  mutate(across(`Compound`, str_replace, "Trimethylamine N-oxide", "TMAO")) %>%
  mutate(across(`Compound`, str_replace, "Betaine", "Glycine betaine"))

fig.1 <- ggplot(dat.fig, aes(x = sample, y = mean.Nmol.C, fill = reorder(Compound, order))) +
  geom_col(alpha = 0.8, color = "black", width = 0.8, size = 0.2) +
  theme_classic() +
  ylab("Concentration (nM C)") +
  xlab("Sample") +
  labs(fill = "Compound") +
  scale_fill_manual(values = lacroix_palette(type = "paired", n = 12)) +
  scale_y_continuous(expand = c(0, NA, NA, 6500), limits = c(0,650)) +
  theme(legend.text=element_text(size=rel(0.5)),
        axis.text=element_text(size=rel(0.5)),
        legend.title = element_text(size=rel(0.7)),
        axis.title=element_text(size=rel(0.7))) +
  theme(legend.position = "bottom")
fig.1

fig.2 <- ggplot(dat.fig, aes(x = sample, y = mean.Nmol.C, fill = reorder(Compound, order))) +
  geom_col(alpha = 0.8, color = "black", width = 0.8, position = "fill", size = 0.2) +
  theme_classic() +
  ylab("Mole Fraction C (%)") +
  xlab("Sample") +
  labs(fill = "Compound") +
  scale_fill_manual(values = lacroix_palette(type = "paired", n = 12)) +
  scale_y_continuous(expand = c(0, NA, NA, NA)) +
  theme(legend.text=element_text(size=rel(0.5)),
        axis.text=element_text(size=rel(0.5)),
        legend.title = element_text(size=rel(0.7)),
        axis.title=element_text(size=rel(0.7))) 
fig.2

leg.1 <- get_legend(fig.1)

ggarrange(ggarrange(fig.1, fig.2, labels = c("A", "B"), legend = "none", nrow = 1, align = "v"), 
          leg.1, nrow = 2, heights = c(4, 1))


ggsave(filename = "Figures/Outputs/Enviro_Conc_CXSPE.pdf", height = 5, width = 5)
ggsave(filename = "Figures/Outputs/Enviro_Conc_CXSPE.jpg", height = 5, width = 5, dpi = 300)






##########################Try looking at 12 most abundant compounds concentration wise
###get 12 most abundant compounds
ma.dat <- dat %>%
  group_by(Compound) %>%
  filter(!Compound == "Salicylic Acid") %>%
  filter(!Compound == "Ornithine") %>%
  summarize(max.abu = max(mean.conc)) %>%
  slice_max(max.abu, n = 11) %>%
  mutate(order = row_number()) %>%
  select(Compound, order)

####Summarize the mean.conc, nmol.C, and nmol.N of the less abundant compounds
dat.other <- full_join(dat, ma.dat) %>%
  filter(is.na(order)) %>%
  group_by(sample) %>%
  summarise(mean.conc = sum(mean.conc),
            mean.Nmol.C = sum(mean.Nmol.C),
            mean.Nmol.N = sum(mean.Nmol.N, na.rm = TRUE)) %>%
  mutate(Compound = "Other",
         order = 12)



###Put most abundant and other data together for figure  
ma.dat.fig <- full_join(dat, ma.dat) %>%
  filter(!is.na(order)) %>%
  select(Compound, sample, mean.conc, mean.Nmol.C, mean.Nmol.N, order)

dat.fig <- rbind(ma.dat.fig, dat.other) 

#Change Trimethylaammonium Propionate (TMAP) to just TMAP
dat.fig <- dat.fig %>%
  ungroup() %>%
  mutate(across(`Compound`, str_replace, "Trimethylammonium Propionate", "TMAP")) %>%
  mutate(across(`Compound`, str_remove, "\\(TMAP\\)")) %>%
  mutate(across(`Compound`, str_replace, "Trimethylamine N-oxide", "TMAO")) %>%
  mutate(across(`Compound`, str_replace, "Betaine", "Glycine betaine"))

fig.3 <- ggplot(dat.fig, aes(x = sample, y = mean.conc, fill = reorder(Compound, order))) +
  geom_col(alpha = 0.8, color = "black", width = 0.8, size = 0.2) +
  theme_classic() +
  ylab("Concentration (nM)") +
  xlab("Sample") +
  labs(fill = "Compound") +
  scale_fill_manual(values = lacroix_palette(type = "paired", n = 12)) +
  scale_y_continuous(expand = c(0, NA, NA, 200), limits = c(0,175)) +
  theme(legend.text=element_text(size=rel(0.5)),
        axis.text=element_text(size=rel(0.5)),
        legend.title = element_text(size=rel(0.7)),
        axis.title=element_text(size=rel(0.7))) +
  theme(legend.position = "bottom")
fig.3

ggsave(filename = "Figures/Outputs/Enviro_Conc_CXSPE_nM_supplemental.pdf", height = 5, width = 5)










