#Exercise 4 
#This is a variation of the code for Exercise 4 plotting results in one figure

#packages
library(deSolve)

#par(mfrow=c(1,2))

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
t <- c(0: 80)            #time vector, h
conc <- read.csv("exposure_scenario.csv", header=FALSE, sep=",", skip=1)

for (i in 1:2){

  
     # if (i==1) conc[,i+1]= conc[,i+1] *3
      
      Cw <- approxfun(t, conc[,i+1],  rule = 2, method = "constant") #interpolation
      
      #Initial state variables
      ini=c(Ci=0)
      
      #parameters
      pars<-c(ku=0.01750034, 
              ke=0.01566554) 
      
      #call function, collect outputs
      out<- lsoda(y = ini,times = t, func = F_tk, pars, input=Cw)
      
      #plot results
      if (i ==1) { 
        # add code here  
           plot(out[,1], out[,2], type="l", lty= 1, lwd=2, xlab = "Time [h]", 
           ylab = "Internal concentration [g/kg]", ylim=c(0, 0.3))
      }  else { 
          # add code here } 
          lines (out[,1], out[,2], type="l", lty= 1, lwd=2, col="red")
      }
      
}
