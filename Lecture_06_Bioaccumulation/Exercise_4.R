#Exercise 4
# Simulate bioaccumulation for a peak exposure scenario and for constant exposure of the area under the curve concentration over time 
# using the parameters estimated in Exercise 3.
# Show both scenarios in a graph.
# What is the difference between peak exposure and constant exposure in terms of bioaccumulation?
  
  
#packages, settings
library(deSolve)
par(mfrow=c(2,2))

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
conc <- read.csv("exposure_scenario.csv", header=FALSE, sep=",", skip=1)
t <- conc[,1]           #time vector

for (i in 1:2){

    Cw <- approxfun(t, conc[,i+1],  rule = 2, method = "constant") #interpolation
    
    #Initial state variables
    ini=c(Ci=0)
    
    #parameters
    pars<-c(ku=0.0175, 
            ke=0.0157) 
    
    #call function, collect outputs
    out<- lsoda(y = ini,times = t, func = F_tk, pars, input=Cw)
    
    #plot results
    plot(conc [,1], conc [,i+1], type="l", lty= 1, lwd=2, xlab = "Time [h]", ylab = "External concentration [g/kg]", ylim=c(0, 2))
    plot(out[,1], out[,2], type="l", lty= 1, lwd=2, xlab = "Time [h]", ylab = "Internal concentration [g/kg]", ylim=c(0, 0.25))

}  