#Exercise 6
# Implement the stochastic death toxicodynamic assumption within the presented TK-TD model 
# by using the scaled damage as dose metric (GUTS RED, slide 17). 
# Therefore, calculate the hazard rate and survival as function of time. Parameter values need to be estimated in the next exercises.
# Try different parameter values for kd, b and z – does the model behave as expected?
# Try different exposure concentrations – does the model behave as expected?
  
  
#packages, settings
library(deSolve)


#toxicokinetic-toxicodynamic function 
F_tktd <- function(t, y,  pars, input) {
  with(as.list(c(y, pars)), {
    #exposure concentration
    Cw     <- input(t)
    #damage
    dD <- kd*(Cw - D)
    #hazard
    dH  <- b*(max(0,D-z))+hb
    
    return(list(c(dD, dH))) 
  })
}

#Input 
conc <- rep(0.3,50)
t <- c(1:length(conc))         #time vector
Cw <- approxfun(t, conc,  rule = 2, method = "constant") #interpolation


#Input from file
# conc <- read.csv("exposure_scenario.csv", header=FALSE, sep=",", skip=1)
# t <- conc[,1]           #time vector
# Cw <- approxfun(t, conc[,2],  rule = 2, method = "constant") #interpolation


#Initial state variables
ini=c(D=0, H=0)

#parameters
pars<-c(z=0.1, 
        b=2,
        kd = 0.0175,
        hb=0) 

#call function, collect outputs
out<- lsoda(y = ini,times = t, func = F_tktd, pars, input=Cw)

#plot results
plot(out[,1], exp(-out[,3]), type="l", lty= 1, lwd=2, col="black", xlab = "Time [h]", ylab = "Survival [-]", ylim=c(0, 1))


    