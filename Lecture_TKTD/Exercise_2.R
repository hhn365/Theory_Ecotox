#Exercise 2
# Bioaccumulation of cadmium has experimentally been measured in Daphnia magna. Elimination [%] has been determined after transfering daphnids to clean media (for data see the xlsx file). 
#  Curate the data file in an R readable format and read the data from file. Present both data (dots) and model simulation (line) in one plot
#   Estimate the parameters kA and kE by trial and error (intially assume kA = kE )

#packages
library(deSolve)

#toxicokinetic function
F_tk <- function(t, y,  pars, input) {
  with(as.list(c(y, pars)), {
    #exposure concentration
    Cw     <- input(t)
    #internal concentration
    dCi <- ku * Cw - ke*Ci 
    return(list(c(dCi))) 
  })
}

#Inputs
t <- c(0: 144)            #time vector
Cw<- rep(0.217, length(t)) #vector of exposure concentrations
Cw[73:144]<-0
Cw <- approxfun(t, Cw,  rule = 2, method = "constant") #interpolation

#Initial state variables
ini=c(Ci=0)

#parameters
 p=0.03
 pars<-c(ku=p,
         ke=p)

 # pars<-c(ku=0.0185, 
 #         ke=0.017) 

#call function, collect outputs
out<- lsoda(y = ini,times = t, func = F_tk, pars, input=Cw)

#plot results
plot(out[,1], out[,2], type="l", lty= 1, lwd=2, xlab = "Time [h]", 
     ylab = "Internal concentration [g/kg]", ylim=c(0, 0.2))

data<- read.csv("bioaccumulation_data_curated.csv", header=TRUE)
 points(data[,1], data[,2])





 #points(data$time,data$concentration)