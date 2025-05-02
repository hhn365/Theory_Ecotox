##------------------------------------------------------------------------------
## GUTS parallel MCMC parameter estimation script (since v2 with Morse package)
## Author: Benoit Goussen (benoit.goussen@gmail.com) (benoit.goussen@ibacon.com)
## Date: 2016/11/30 (initial)- 2018/03/01
## Version: 2.1
## Packages used: Morse, tidyr, ggplots2, JAGS
## Copyright 2016, 2017, 2018 Benoit Goussen
## Licence: GPL-3.0+
##------------------------------------------------------------------------------
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

##------------------------------------------------------------------------------

## Modified by Andre Gergs for teaching purposes


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ INITIALISATION ~~~~ -----------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## * Clear all and set seed  ===================================================
rm(list = ls()) # Clear all objects
graphics.off() # Close all graphics
set.seed(42) # Set the seed for the random number generators

# * Check if packages are installed and load ==================================
 if (!require(tidyr)) {
   install.packages("tidyr")
   library(tidyr)
 }
 if (!require(ggplot2)) {
   install.packages("ggplot2")
   library(ggplot2)
 }
 if (!require(morse)) {
   install.packages("morse")
   library(morse)
 }

packageVersion('morse')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ DATA ~~~~ ---------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## * Load and work on data =====================================================
## load RData if calibration was already performed
load("resultsTest-daphnia.RData")

### ** Read in data ####
dfDat <- read.csv("dataTest_daphnids_Cd.csv", header = T, sep = ",", dec = ".") # Training data

# Example format of the dataset. Data at time t=0 must be provided
# Id  time  conc  Nsurv replicate
# 1   0     0     7     1
# 2   2     0     7     1
# 3  24     0     7     1
# 4  48     0     7     1
# 5  72     0     7     1
# 6  96     0     7     1
# 7   0     0.5   7     2
# 8   2     0.5   7     2
# 9  24     0.5   7     2
# 10 48     0.5   7     2
# 11 72     0.5   7     2
# 12 96     0.5   7     2
# 13  0     1.0   7     3
# 14  2     1.0   7     3
# ...


## * Create survData object ====================================================
# Check that the data are in the right format and create the survData object
survDataCheck(dfDat)
dat <- survData(dfDat)


## * Plot the data and a summary ===============================================
plot(dat, pool.replicate = TRUE, ylab = "Number of mobile individuals")
#plotDoseResponse(dat, addlegend = TRUE, log = T, ylab = "Mobility rate")
summary(dat)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ INFERENCE AND PLOTS ~~~~ ------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## * FIt the GUTS SD and make plots ============================================
fit_SD <- survFit(dat, quiet = FALSE, model_type = "SD", dic.compute = T) 
summary(fit_SD) # Summary of the inference

# Plot PPC
gg1 <- ppc(fit_SD)
gg1 + geom_abline(intercept = 0, slope = 1, linetype = "dashed")

# Perform prediction using the time and concentrations of the evaluation dataset
 gg <- plot(fit_SD, xlab = "Time [d]", ylab = "Survival [-]", adddate = TRUE)# Plot inference results
 gg + geom_point(data = dat,
                 aes(x = time, y = Nsurv / Ninit, group = replicate),
                 col = "black", pch = 20, size = 3
 )+ theme(strip.background= element_rect(fill = "white")) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ Performance Critera ~~~~ ------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
predict_SD <- predict_Nsurv(fit_SD, data_predict = dat)
predict_Nsurv_check(predict_SD)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ FINAL ~~~~ --------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
save.image("resultsTest-daphnia.RData") # Save workspace
write.table(summary(fit_SD),sep=',',dec='.',"GUTS_SD_parameters.csv",row.names=FALSE , col.names=TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~ END ~~~~ ----------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
