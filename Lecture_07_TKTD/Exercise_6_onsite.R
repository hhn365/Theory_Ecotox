#Exercise 5
# Replace the bioaccumulation model by the one-parameter function for the scaled internal concentration Ci.
# Compare modeling results for the scaled internal concentration with the ones obtained from Exercise 4
  
  
#packages, settings
library(deSolve)
par(mfrow=c(1,2))

#toxicokinetic function
F_tk <- function(t, y,  pars, input) {
  with(as.list(c(y, pars)), {
    #exposure concentration
    Cw     <- input(t)
    #damage
    dD <- kd*(Cw - D)
    #hazard
    dH <- b*max(D-z,0)+hb
    
    return(list(c(dD,dH))) 
  })
}


#Inputs
conc <- read.csv("exposure_scenario.csv", header=FALSE, sep=",", skip=1)
t <- conc[,1]           #time vector

#for (i in 1:2){
i=1

    Cw <- approxfun(t, conc[,i+1],  rule = 2, method = "constant") #interpolation
    
    #Initial state variables
    ini=c(D=0, H=0)
    
    #parameters
    pars<-c(kd =0.001,
            b=10,
            z=0.01,
            hb=0
            
            ) # 0.0166
    
    #call function, collect outputs
    out<- lsoda(y = ini,times = t, func = F_tk, pars, input=Cw)
    
    
    #plot results
    plot(out[,1], out[,2], type="l", lty= 1, lwd=2, col="orange", xlab = "Time [h]", ylab = "Damage [g/kg]", ylim=c(0, 0.25))
    plot(out[,1], exp(- out[,3]), type="l", lty= 1, lwd=2, col="red", xlab = "Time [h]", ylab = "survival [g/kg]", ylim=c(0, 1))
    
#}  