#
#
#
#
#
#
#
#
library(tidyverse)
library(readr)
library(ggpubr)


###Define inputs
CXC_PPL_data_file <- "Intermediates/Analytical_Validation/AV_CX_PPL_HQ_combined.csv"
Std.file <- "Meta_Data/Ingalls_Standards/Ingalls_Lab_Standards_NEW.txt"

#Load in data:
All.dat <- read_csv(CXC_PPL_data_file)

###Determine Overlaps 
####Find overlaps 
overlap.CX <- All.dat %>%
  select(MF, method, Fraction) %>%
  unique() %>%
  filter(method == "CX-SPE") %>%
  mutate(score.cat = 1) %>%
  select(MF, Fraction, score.cat)

overlap.PPL <- All.dat%>%
  select(MF, method, Fraction) %>%
  unique() %>%
  filter(method == "PPL-SPE") %>%
  mutate(score.ppl = 2) %>%
  select(MF, Fraction, score.ppl)

overlap.dat <- full_join(overlap.CX, overlap.PPL) %>%
  rowwise() %>%
  mutate(overlap.score = sum(score.cat, score.ppl, na.rm = TRUE)) %>%
  mutate("Overlap" = case_when(overlap.score == 1 ~ "CX_only", 
                               overlap.score == 2 ~ "PPL_only",
                               overlap.score == 3 ~ "Both"
  )) %>%
  select(MF, Fraction, Overlap)

####
#
#load in standards
std.info <- read_csv(Std.file)

std.form <- std.info$Emperical.Formula

formula2elements <- function(formula_vec){
  split_formulas <-   regmatches(gregexpr(formula_vec, pattern = "[A-Z][a-z]*[0-9]*"), x = formula_vec)
  elements_in <- lapply(split_formulas, gsub, pattern = "[0-9]", replacement = "")
  element_counts <- lapply(split_formulas, gsub, pattern = "[A-z]", replacement = "")
  element_counts <- lapply(element_counts, function(x){
    x[!nchar(x)]<-1
    return(as.numeric(x))
  })
  mapply(`names<-`, element_counts, elements_in, SIMPLIFY = FALSE)
}

form_table <- formula2elements("C6H12O6")[[1]]
form_df <- dplyr::bind_rows(formula2elements(c("C6H12O6", "C5H11NO")))

std.form.df <- dplyr::bind_rows(formula2elements(std.form)) %>%
  select(C, H, O)

std.elements <- bind_cols(std.info, std.form.df) %>%
  select(Compound.Name_old, Emperical.Formula, C, H, O) %>%
  rename(MF = Compound.Name_old)



#empir.std <- std.info %>%
#  filter(Priority = TRUE) %>%
#  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
 #                   str_extract(Emperical.Formula, "^C\\d"), 
#                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
#  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
#  
#  mutate(H = ifelse(is.na(str_extract(Emperical.Formula, "^H\\d\\d")),
#                    str_extract(Emperical.Formula, "^H\\d"), 
#                    str_extract(Emperical.Formula, "^H\\d\\d"))) %>%
#  mutate(N = ifelse(str_detect(Emperical.Formula, "^N\\D"),
#                    1, str_extract(Emperical.Formula, "^N\\d")))%>%
#  select(Compound.Name_old, Emperical.Formula, C, H, N)
  
  
  
  
  
  
  
  
 # mutate(H = case_when(str_detect(.$Emperical.Formula, "H\\d\\d") ~ str_extract(.$Emperical.Formula, "^H\\d\\d"),
 #                      str_detect(.$Emperical.Formula, "^H\\d\\D") ~ str_extract(.$Emperical.Formula, "^H\\d"))) %>% 
#  mutate(C = case_when(str_detect(.$Emperical.Formula, "^C\\d\\d") ~ str_extract(.$Emperical.Formula, "^C\\d\\d"),
 #                      str_detect(.$Emperical.Formula, "^C\\d") ~ str_extract(.$Emperical.Formula, "^C\\d"))) %>% #,
#                       #str_detect(.$Emperical.Formula, "^C\\D") ~ 1)) %>%
 # mutate(H = case_when(str_detect(.$Emperical.Formula, "^H\\d\\d") ~ str_extract(.$Emperical.Formula, "^H\\d\\d"),
#                       str_detect(.$Emperical.Formula, "^H\\d\\D") ~ str_extract(.$Emperical.Formula, "^H\\d"))) %>% 
#  select(Compound.Name_old, Emperical.Formula, C, H)

  
  
  
  
#  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
 #                   str_extract(Emperical.Formula, "^C\\d"), 
#                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
#  mutate(H = ifelse(is.na(str_extract(Emperical.Formula, "^H\\d\\d")),
#                    str_extract(Emperical.Formula, "^H\\d"))) %>%
#  mutate(O = ifelse(str_detect(Emperical.Formula, "O\\D"),
 #                   1, str_extract(Emperical.Formula, "O\\d"))) %>%
##  mutate(O = as.numeric(str_replace_all(O, "O", ""))) %>%
#  select(Compound.Name_old, Emperical.Formula, C, H, O)
#

##parse #s of C, H, and O from empirical formula
empir.std <- std.info %>%
  filter(Priority = TRUE) %>%
  mutate(C = ifelse(is.na(str_extract(Emperical.Formula, "^C\\d\\d")),
                    str_extract(Emperical.Formula, "^C\\d"), 
                    str_extract(Emperical.Formula, "^C\\d\\d"))) %>%
  mutate(C = as.numeric(str_replace_all(C, "C", ""))) %>%
  mutate(N = ifelse(str_detect(Emperical.Formula, "N\\D"),
                    1, str_extract(Emperical.Formula, "N\\d")))%>%
  mutate(N = as.numeric(str_replace_all(N, "N", ""))) %>%
  mutate(O = ifelse(str_detect(Emperical.Formula, "O\\D"),
                    1, str_extract(Emperical.Formula, "O\\d"))) %>%
  mutate(O = as.numeric(str_replace_all(O, "O", ""))) %>%
  mutate(S = ifelse(str_detect(Emperical.Formula, "S\\D"),
                    1, str_extract(Emperical.Formula, "S\\d"))) %>%
  mutate(S = as.numeric(str_replace_all(S, "S", ""))) %>%
  mutate(P = ifelse(str_detect(Emperical.Formula, "P\\D"),
                    1, str_extract(Emperical.Formula, "P\\d"))) %>%
  mutate(P = as.numeric(str_replace_all(P, "P", ""))) %>%
  mutate(H = ifelse(str_detect(Emperical.Formula, "H\\d\\d"),
                    str_extract(Emperical.Formula, "H\\d"), 1)) %>%
  mutate(H = as.numeric(str_replace_all(H, "H", ""))) %>%
  select(Compound.Name_old, Emperical.Formula, C, H, O, N, P, S, m.z, RT..min.) %>%
  rename(MF = Compound.Name_old)

#calculate C/O, C/H, and C/N ratios
empir.dat <- left_join(overlap.dat, std.elements) %>%
  unique() %>%
  mutate(`O/C` = O/C,
         `H/C` = H/C) 


library(wesanderson)
library(viridis)

###VMake Van Krevelen Plot
VK.plot <- ggplot(empir.dat, aes(x = `O/C`, y=`H/C`, color = Overlap)) +
  geom_jitter(size = 2, alpha = 0.5, height = 0.0075, width = 0.0075) +
  theme_classic() + 
 # scale_color_brewer(palette = "Dark2") +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  labs(color = "Compound Classification") +
  theme(axis.text= element_text(size = 6))

#Make C/H Density Plot
CH.den.plot <- ggplot(empir.dat, aes(x = `H/C`, fill = Overlap, color = Overlap)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  facet_grid(rows = vars(Overlap)) +
  theme(strip.text = element_text(size = 8),
        axis.text= element_text(size = 6),
        axis.line = element_line(size = 0.5),
        strip.background = element_rect(size = 0, color = "white"),
        legend.position = "none")

#Make C/O density plot
CO.den.plot <- ggplot(empir.dat, aes(x = `O/C`, fill = Overlap, color = Overlap)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  scale_fill_manual(values = wes_palette("Darjeeling1")) +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  facet_grid(rows = vars(Overlap)) +
  theme(strip.text = element_text(size = 8),
        axis.text= element_text(size = 6),
        axis.line = element_line(size = 0.5),
        strip.background = element_rect(size = 0, color = "white"),
        legend.position = "none")

###
ggarrange(VK.plot, labels = c("A"),
          ggarrange(CH.den.plot, CO.den.plot, labels = c("B", "C"), ncol = 2), nrow = 2, common.legend = TRUE, align = "v")
#
ggsave(filename = "Figures/Outputs/VanKrevelen_plot.pdf", height = 4.5, width = 5)


















