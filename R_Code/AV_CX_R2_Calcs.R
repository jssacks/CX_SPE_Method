#
#
#

#
#This Script calculates R2 values from the calibration curve made in ALOHA seawater for all compounds to evaluate the linearity
# of the method response to nM level changes in analyte concentration


#packages
library(tidyverse)
library(readr)

#source
source("R_Code/CXSPE_Functions.R")

##Define Inputs
CX.file <- "Intermediates/Analytical_Validation/AV_CX_Full_Dat.csv"


###
CX.dat <- read_csv(CX.file)

###Select data for calibration curve
cal.dat <- CX.dat %>%
  filter(samp == "A" |
           samp == "LinA") %>%
  mutate(Spike.Conc = case_when(
    replicate == "0point1" ~ 0.1,
    replicate == "0point5" ~ 0.5,
    replicate == "SBe" ~ 1))

######Create calibration curve to determine R2 values

###make data frame of just MFs
dat.MF <- cal.dat %>%
  ungroup()%>%
  select(MF) %>%
  unique()

####Run Calibration curve function to get R2 output for every compound
for (i in seq_along(dat.MF$MF)) {
  res <- dat.MF$MF[i]
  cal.calc <- Calibcurv_calc(cal.dat, res)
  res.df <- data.frame(cal.calc)
  out <- full_join(out, res.df)
}
cal.curv.output <- out

###write data to intermediates file
write_csv(cal.curv.output, path = "Intermediates/Analytical_Validation/AV_CX_R2_Dat.csv")



