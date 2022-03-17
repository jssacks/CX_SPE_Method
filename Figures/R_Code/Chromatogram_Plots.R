#####
#####
## This script creates the extracted ion chromatogram (EIC) plots showing A) different concentrations 
# of compounds extracted using CX-SPE (main text figure 1) and B) the removal of contaminant glycine betaine
# using the CXC approach (supplementary figure 1)

###Load in packages
library(tidyverse)
library(RaMS)
library(data.table)

# Main Text Figure 1 ------------------------------------------------------

#define .mzml file inputs
ms.blk <- "Raw_Data/MZML_Files/210720_Blk_Blank1_B.mzML"
ms.ns <- "Raw_Data/MZML_Files/210720_Smp_DS_NoSPike_A.mzML"
ms.1nm <- "Raw_Data/MZML_Files/210720_Smp_DS_1nM_C.mzML"
ms.5nm <- "Raw_Data/MZML_Files/210720_Smp_DS_5nM_C.mzML"
ms.25nm <- "Raw_Data/MZML_Files/210720_Smp_DS_25nM_C.mzML"

###load in mass spec files
ms.dat <- grabMSdata(c(ms.blk, ms.ns, ms.1nm, ms.5nm, ms.25nm))

###Extract out homarine, GBT, and TMAO
homarine.mass <- 138.05503
GBT.mass <- 118.086804
TMAO.mass <- 76.076
DMSAc.mass <- 121.032327

mzs.eic <- c(`Homarine and Trigonelline`=homarine.mass, `Gylcine Betaine and Valine`=GBT.mass, TMAO=TMAO.mass, `DMS-Ac`=DMSAc.mass)
eic.dat <- imap_dfr(mzs.eic, function(mz_i, name){
  cbind(ms.dat$MS1[mz%between%pmppm(mz_i, ppm=10)], name)
})

###Clean up data file and remove unneaded data outside RT window
eic.dat.2 <- eic.dat %>%
  filter(!rt>10) %>%
  filter(!rt<4.5) %>%
  mutate(Spike_Concentration = case_when(str_detect(.$filename, "Blk") ~ "Blank",
                                         str_detect(.$filename, "NoSPike") ~ "No Spike",
                                         str_detect(.$filename, "_1nM") ~ "1nM",
                                         str_detect(.$filename, "_5nM") ~ "5nM",
                                         str_detect(.$filename, "_25nM") ~ "25nM",)) %>%
  select(-filename)

###Add additional DMS-Ac data as there is no signal at this mass for 
# most of the chromatograms (especially the Blank)
rt <- c(450:804, 950:1000)
DMSAc.dat <- data.frame(rt) %>%
  mutate(rt = rt/100,
         name = "DMS-Ac",
         int = 5000)

Spike_Concentration <- c("Blank", "No Spike", "1nM", "5nM", "25nM")
name <- c("DMS-Ac", "DMS-Ac", "DMS-Ac", "DMS-Ac", "DMS-Ac")
DMSAc.filenames <- tibble(Spike_Concentration, name)

DMSAc.Add <- full_join(DMSAc.filenames, DMSAc.dat)  

##Add supplemental DMS-Ac data into the overall data frame
eic.dat.3 <- rbind(eic.dat.2, DMSAc.Add, fill = TRUE) %>%
  mutate(Spike_Concentration = as.factor(Spike_Concentration))

##Add in ordering to Spike Concentration
eic.dat.3$Spike_Concentration <- ordered(eic.dat.3$Spike_Concentration, 
                                         levels = c("25nM",
                                                    "5nM",
                                                    "1nM",
                                                    "No Spike",
                                                    "Blank"))
  
###Define Palette 
cbPalette <- c("#0072B2", "#CC79A7", "#009E73", "#E69F00", "#999999")
lab.df <- tibble(x = c(5.5, 9, 6, 8.25),
                 y = c(1.5e6, 3e6, 4e7, 1e7),
                 text = c("Homarine", "Trigonelline", "Glycine Betaine", "Valine"),
                 name = c("Homarine and Trigonelline", "Homarine and Trigonelline", "Gylcine Betaine and Valine", "Gylcine Betaine and Valine"))

###Create Plot and save as a pdf
Chrom.Plot.1 <- ggplot(eic.dat.3)  +
  geom_smooth(aes(x=rt, y=int, color = Spike_Concentration), size = 0.9, se=FALSE, span = 0.05, linetype = 1) +
  facet_wrap(~name, ncol = 2, scales = "free_y") +
  theme_classic() +
  labs(color = "Spike Concentration") + 
  scale_color_manual(values = cbPalette) +
  xlab("Retention Time (min)") +
  ylab("Intensity") +
  geom_text(data = lab.df, aes(x = x, y = y, label = text)) +
  theme(legend.position = c(.15, .8), 
        legend.text = element_text(size = 12),
        legend.background = element_rect(linetype = 1, size = 0.25, color = 1),
        strip.background = element_rect(color = "white"),
        strip.text.x = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        panel.spacing.y = unit(1.5, "lines"))
Chrom.Plot.1
ggsave("Figures/Outputs/MainTextFig1_chomatogram.pdf")



# Supplementary Figure 1 - CXC Chromatograms -------------------------------

###load in mass spec files
ms.MQ <- "Raw_Data/mzML_Files/210209_Smp_SPE_Wash_MQ_A.mzML"
ms.W1 <- "Raw_Data/mzML_Files/210209_Smp_SPE_Wash_W1_A.mzML"

###Get MS Data
ms.dat.wash <- grabMSdata(c(ms.MQ, ms.W1))

###Clean and organizae MS Data, grab just data for glycine betaine (GBT)
wash.gbt <- ms.dat.wash$MS1[mz%between%pmppm(GBT.mass, ppm = 10)] %>%
  mutate(Sample = case_when(str_detect(.$filename, "MQ") ~ "MQ-H2O", 
                            str_detect(.$filename, "W1") ~ "CXC-H2O"))

###Make Figure and Save to Outputs
CXC.fig <- ggplot(wash.gbt) +
  geom_line(aes(x=rt, y=int, color=Sample), size=1, alpha = 0.8) +
  theme_classic() +
  xlab("Retention Time (min)") +
  ylab("Intensity")
ggsave("Figures/Outputs/SupFig1_CXC_chromatogram.pdf")
