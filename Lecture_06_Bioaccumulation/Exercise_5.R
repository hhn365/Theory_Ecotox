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
    #internal concentration
    dCi <- ku * Cw - ke*Ci 
    dCi2 <- kd*(Cw - Ci2)
    return(list(c(dCi, dCi2))) 
  })
}


#Inputs
conc <- read.csv("exposure_scenario.csv", header=FALSE, sep=",", skip=1)
t <- conc[,1]           #time vector

for (i in 1:2){

    Cw <- approxfun(t, conc[,i+1],  rule = 2, method = "constant") #interpolation
    
    #Initial state variables
    ini=c(Ci=0, Ci2=0)
    
    #parameters
    pars<-c(ku=0.0175, 
            ke=0.0157,
            kd =0.0175  ) # 0.0166
    
    #call function, collect outputs
    out<- lsoda(y = ini,times = t, func = F_tk, pars, input=Cw)
    
    
    #plot results
    plot(out[,1], out[,2], type="l", lty= 1, lwd=2, col="orange", xlab = "Time [h]", ylab = "Internal concentration [g/kg]", ylim=c(0, 0.25))
     lines(out[,1], out[,3], lty= 1, lwd=2, col="steelblue")
       legend(1,0.25, lty = 1,lwd=2, col= c("steelblue","orange"), legend=c("1 par", "2 pars"), box.col = "white") 

}  