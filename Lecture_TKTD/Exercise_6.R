#Exercise 5
# Replace the bioaccumulation model by the one-parameter function for the scaled internal concentration Ci.
# Compare modeling results for the scaled internal concentration with the ones obtained from Exercise 4
  
  
#packages, settings
library(deSolve)
#par(mfrow=c(1,2))

#toxicokinetic function
F_tktd <- function(t, y,  pars, input) {
  with(as.list(c(y, pars)), {
    #exposure concentration
    Cw     <- input(t)
    #internal concentration
     
    dD <- kd*(Cw - D)
    dH <- b * max(D-z,0)+hb
    return(list(c(dD, dH))) 
  })
}


#Inputs
conc <- read.csv("exposure_scenario.csv", header=FALSE, sep=",", skip=1)
t <- conc[,1]           #time vector

for (i in 1:1){

    Cw <- approxfun(t, conc[,i+1],  rule = 2, method = "constant") #interpolation
    
    #Initial state variables
    ini=c(D=0, H=0)
    
    #parameters
    pars<-c(z=0.01, #threshold
            b=0.5,  #slope
            hb=0,   #background hazard rate
            kd =0.0175 ) # dominant rate
    
    #call function, collect outputs
    out<- lsoda(y = ini,times = t, func = F_tktd, pars, input=Cw)
    
    
    #plot results
    plot(out[,1], exp(-out[,3]), type="l", lty= 1, lwd=2, col="orange", xlab = "Time [h]", ylab = "Internal concentration [g/kg]", ylim=c(0, 1))
     #lines(out[,1], out[,3], lty= 1, lwd=2, col="steelblue")

}  