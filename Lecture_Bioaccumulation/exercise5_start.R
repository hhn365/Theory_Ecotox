library(deSolve)

#implement and call a function in R (‚deSolve‘ package)
my_function<- function(t, y,  pars, input) {
  with(as.list(c(y, pars)), { 
    Cw <- input(t)
    #change in state
    dCi=ku*Cw-ke*Ci  
    dCi2=kd*(Cw-Ci2)
    return(list(c(dCi, dCi2)))   })}

#input
t <- c(0: 80)            #time vector, h
conc <- read.csv("exposure_scenario.csv", header=FALSE, sep=",", skip=1)

Cw <- approxfun(t, conc[,2],  rule = 2, method = "constant") #interpolation

#initial state variables
yini<- c(Ci=0, Ci2=0)
#parameters
pars<-c(
  ku=0.01750034,
  ke=0.01566554,
  kd = 0.01750034
)

#call function, collect outputs
output<- lsoda(y = yini, times = t, func = my_function, pars, input=Cw)
                   

#plot results
plot(output[,1], output[,2], type="l", lty= 1, lwd=2, xlab = "Time [h]", 
     ylab = "Internal concentration [g/kg]")
       lines(output[,1], output[,3], lty= 1, lwd=2, col="steelblue")

       legend(70,32 , lty = 1,lwd=2, col= c("steelblue","black"), legend=c("high", "low"), box.col = "white") 