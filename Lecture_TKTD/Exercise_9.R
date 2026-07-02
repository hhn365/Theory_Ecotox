#Exercise 9
# Calculate LP50 values for the two exposure profiles.
# Do you observe any differences between the exposure profiles? Why?

  
#packages, settings
library(deSolve)
library(drc)



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

#parameters
pars<-c(z=7.16E-02, 
        b=1.70E-02, #kk in morse
        kd=5.23E-03,
        hb=0) 

#results matrix
no_factors=100
res<-matrix(NA, nrow = no_factors, ncol=3)

for (i in 1:2) {
    
  for (f in 1:no_factors) {
      Cw <- approxfun(t, conc[,1+i]*f,  rule = 2, method = "constant") #interpolation
      
      #Initial state variables
      ini=c(D=0, H=0)

      #call function, collect outputs
      out<- lsoda(y = ini,times = t, func = F_tktd, pars, input=Cw)
     
     res[f,1] <- f #loop counter is used as multiplication factor 
     res[f,1+i] <- 1- exp(-out[length(t),3]) #effect (1-survival) at the end saved
  }
   
} 

#data frame
res<-as.data.frame(res)
colnames(res)<-c("factor", "constant", "variable")

par(mfrow=c(1,2))

#fit dose response models for Lpx:  
#constant
model1 <- drm(res [,2]~res [,1], fct=LL.2(), type = "binomial")
 plot(model1,  main='Constant exposure',log='x',  ylim=c(0,1),ylab = "Effect [-]",xlab = "Multiplication factor [#]", pch=19,cex = 1.3)

#variable
model2 <- drm(res [,3]~res [,1], fct=LL.2(), type = "binomial")
 plot(model2,  main='Variable exposure',log='x',  ylim=c(0,1),ylab = "Effect [-]",xlab = "Multiplication factor [#]", pch=19,cex = 1.3)

#results output
summary(model1)
summary(model2)
