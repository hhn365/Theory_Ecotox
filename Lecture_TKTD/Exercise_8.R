#Exercise 8
# Based on the script from Exercise 6 and the parameter values estimated in Exercise 7 (but hb = 0), 
# simulate the effects for the two exposure profiles used in Exercise 2. 
# What do you observe, and why?
  
  
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

#Input from file
conc <- read.csv("exposure_scenario.csv", header=FALSE, sep=",", skip=1)
t <- conc[,1]           #time vector

par(mfrow=c(2,2))

for (i in 1:2) {
    Cw <- approxfun(t, conc[,1+i],  rule = 2, method = "constant") #interpolation
    
    #Initial state variables
    ini=c(D=0, H=0)
    
    #parameters
    pars<-c(z=0.0716, 
            b=0.017, #kk in morse
            kd=0.00523,
            hb=0) 
    
    #call function, collect outputs
    out<- lsoda(y = ini,times = t, func = F_tktd, pars, input=Cw)
    
    #plot results
    plot(out[,1], out[,2], type="l", lty= 1, lwd=2, col="black", xlab = "Time [h]", ylab = "Scaled damage [mg/L]", ylim=c(0, 0.2))
    
    plot(out[,1], exp(-out[,3]), type="l", lty= 1, lwd=2, col="black", xlab = "Time [h]", ylab = "Survival [-]", ylim=c(0, 1))
    
}    
        