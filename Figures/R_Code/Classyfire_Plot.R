#
#
#
#
#
####
library(tidyverse)
library(readr)
library(ggpubr)


####Define Inputs
Classy_file <- "Meta_Data/Ingalls_Standards/Classy_Ingalls_Stds.csv"
Std_file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"
CXC_PPL_data_file <- "Intermediates/Analytical_Validation/AV_CX_PPL_HQ_combined.csv"
Classy_addition_file <- "Meta_Data/Ingalls_Standards/Classyfire_additional_compounds.csv"
###Load and Combine
Classy.dat <- read_csv(Classy_file) %>% 
  unique() %>%
  select(!CHEMONT) %>%
  pivot_wider(id_cols = key, names_from = Level, values_from = Classification) %>%
  rename(InChI.Key.Name = key)
Std.dat <- read_csv(Std_file) %>%
  mutate(MF = Compound.Name_old,
         column = Column,
         RT = RT..min.)

###
CXC.PPL.dat <- read_csv(CXC_PPL_data_file) %>%
  mutate(good = "Yes") %>%
  mutate(z = ifelse(Fraction == "Neg", -1, 1)) %>%
  mutate(column = ifelse(Fraction == "RP", "RP", "HILIC"))


####
class.add.dat <- read_csv(Classy_addition_file) %>%
  select(-`...6`) %>%
  rename(MF = Compound,
         kingdom = Kingdom,
         superclass = Superclass,
         class = Class,
         subclass = Subclass) %>%
  mutate(MF = str_replace(.$MF, "(3-Carboxypropyl)trimethylammonium (TMAB)",	"(3-Carboxypropyl)trimethylammonium (TMAB)"))


####Combine data
std.class.dat <- full_join(Std.dat, Classy.dat, by = "InChI.Key.Name") %>%
  filter(Priority == "TRUE") %>%
  select(MF, column, z, kingdom, superclass, class, subclass, `level 5`, `level 6`, `level 7`, `level 8`)
full.dat <- full_join(CXC.PPL.dat, std.class.dat, by = c("MF", "z", "column"))


###Just good data
good.dat <- full.dat %>%
  filter(good == "Yes") %>%
  mutate(good = if_else(good == "Yes", "yes", "no", missing = "no")) %>%
  select(good, MF, column, z, kingdom, superclass, class, 
         subclass, `level 5`, `level 6`, `level 7`, `level 8`, Overall.Mean.EE) %>%
  unique() %>%
  filter(!str_detect(MF, ", ")) %>%
  mutate(z = as.factor(z)) %>%
  filter(!z == "2") %>%
  filter(!column == "LipidC18")

###Get good compounds
CXC.PPL.classy.dat <- left_join(CXC.PPL.dat, full.dat) %>%
  select(good, MF, method, column, z, kingdom, superclass, class, 
         subclass, `level 5`, `level 6`, `level 7`, `level 8`, Overall.Mean.EE) %>%
  unique() %>%
  filter(Overall.Mean.EE <= 150)

###Add in compounds that had to be manually annotated 
### TMAB removed because of bad/nonsensical ID by classyfired
additions <- left_join(class.add.dat, CXC.PPL.classy.dat, by = "MF") %>%
  rename(kingdom = kingdom.x,
         superclass = superclass.x,
         class = class.x,
         subclass = subclass.x) %>%
  select(good, MF, method, column, z, kingdom, superclass, class, 
         subclass, `level 5`, `level 6`, `level 7`, `level 8`, Overall.Mean.EE) %>%
  unique() %>%
  filter(!str_detect(.$MF, "TMAB"))

###Add additions into larger dataframe
CXC.PPL.classy.dat.2 <- CXC.PPL.classy.dat %>%
  filter(!MF %in% additions$MF)

CXC.PPL.classy.dat.3 <- rbind(CXC.PPL.classy.dat.2, additions) %>%
  filter(!str_detect(.$MF, "TMAB")) 
  


####Make Plots:

#Superclass Plot_________________

#remove NAs
p1.dat <- CXC.PPL.classy.dat.3 %>%
  filter(!is.na(.$superclass))

plot.1 <- ggplot(p1.dat, aes(y = fct_rev(fct_infreq(superclass)), fill = method)) +
  geom_bar(position = "dodge", width = 0.7, alpha = 0.85) + 
  theme_classic() + 
  scale_fill_manual(values = c("darkred", "darkblue")) +
  scale_x_continuous(expand = c(0, NA)) +
  ylab("Superclass")


#Class Plot_______________________

#remove NAs
p2.dat <- CXC.PPL.classy.dat.3 %>%
  filter(!is.na(.$class))

plot.2 <- ggplot(p2.dat, aes(y = fct_rev(fct_infreq(class)), fill = method)) +
  geom_bar(position = "dodge", width = 0.7, alpha = 0.85) + 
  theme_classic() + 
  scale_fill_manual(values = c("darkred", "darkblue")) +
  scale_x_continuous(expand = c(0, NA)) +
  ylab("Class")

####
ggarrange(plot.1, plot.2, labels = c("A", "B"), common.legend = TRUE, nrow = 2, heights = c(2,5), widths = 5, align = "v")
#save
ggsave(filename = "Figures/Outputs/Classyfire.pdf", height = 6, width = 5)


