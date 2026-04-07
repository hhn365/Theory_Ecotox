#Exercise 1
# Implement the one-compartment bioaccumulation model (equation on slide 18) in R and present model predictions in a graph (line plot)
# Try different parameter values for ku and ke – does the model behave as expected?
# Try different exposure concentrations (including an uptake phase and a depuration phase) – does the model behave as expected?
  
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
t <- c(0: 48)            #time vector
Cw<- rep(0.2, length(t)) #vector of exposure concentrations
Cw[24:48]<-0
Cw <- approxfun(t, Cw,  rule = 2, method = "constant") #interpolation

#Initial state variables
ini=c(Ci=0)

#parameters
pars<-c(ku=0.5, 
        ke=0.4) 

#call function, collect outputs
out<- lsoda(y = ini,times = t, func = F_tk, pars, input=Cw)

#plot results
plot(out[,1], out[,2], type="l", lty= 1, lwd=2, xlab = "Time [h]", 
     ylab = "Internal concentration [g/kg]", ylim=c(0, 0.3))
