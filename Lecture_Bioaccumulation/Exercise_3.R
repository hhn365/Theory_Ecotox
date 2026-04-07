#Exercise 3
# Use the cadmium bioaccumulation data for parameter estimation
#  Estimate the parameters kA and kE by fitting the one-compartment bioaccumulation model to the data using the rbioacc package in R
#  First, curate the data in accordance with the package requirements. For an example see Table 1 in Ratier et al. 2022
#  Read data in R and execute the fitting routine of the package (you may find the code example on the next slide hepful)
#  Plot fitting results and get parameter estimates (see code example)
#  Determine the bioconcentration factor from kinetics (BCF =kA/kE)


#install.packages("rbioacc")
library("rbioacc")

data <- read.csv("bioaccumulation_data_bioacc.csv", header = T, sep = ",", dec = ".") # data

### specify the time at which the accumulation phase ends (here 72 hours)
modeldata <- modelData(data, time_accumulation = c(72))

### fit the TK model
m <- fitTK(modeldata)

### plot the fitting results
plot(m)
q<-quantile_table(m)
q[1:2,]


#calculate BCF----------------------------------------------------------
BCFk_all <- bioacc_metric(m, "k") # "k" for kinetic BCF
### get a summary of the kinetic BCF
for(i in 1:ncol(BCFk_all)){
  BCFk <- quantile(BCFk_all[, i], c(0.5, 0.025, 0.975))
}

### plot the BCF posterior probability distribution
plot(BCFk_all)
### get the TK parameter estimates



