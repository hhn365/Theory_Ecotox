# Day 4
# Implementation of an simplified IBM based on the DEBtox model (Jager and Zimmer 2012) to simulate population dynamcis under chemical exposure.
# We assume variable food and competition, but constant temperature 
# Implementation for teaching purposes only.

#simulation settings--------------------------------------------------------------------------------
tmax    <-   70                   # simulation time [days]
MC.no   <-   10

# number of Monte-Carlo executions
res <- matrix(data = NA, nrow = tmax, ncol = 5) # results matrix
res.pop <- matrix(data = NA, nrow = tmax, ncol = MC.no) # results matrix

#total food availability mg/d
food.daily <-0.5
pop.feed<- rep(c(NA),tmax)

#  Volume of environment in mL
env.size         <- 1000  

# exposure concentration 
Cw      <-rep(c(0), tmax )                      
#Cw [30:35]<-0.16
#Cw [5:9]<-0.2

pMoA <- 1 #Enter physiological mode of action considered for the chemCal effect, here are the ones that are implemented: 
         # 1: feeding/assimilation reduced 
         # 2: Maintanance costs increased 
         # 3: Cost for structure and maturation increased
         # 4: Hazard during oogenesis/embryogenesis

pop.ini <-5   #initial population size (total)
pop.adults<-0 #number of adults in initial population 


# Parameters-----------------------------------------------------------------------------------
v       <- 2.2    # energy conductance
g       <- 0.92   # energy investment ratio          
km      <- 0.51   # somatic maintencance rate coefficalent
Lb      <- 0.86   # length at birth
Lp      <- 2.8    # length at puberty
R_m     <- 37.8   # maximum reproduction rate


s_shrink<- 0.8    # shrinking coefficient for hazard due to starvation
Jxm     <- 0.02503 # scaled maximum ingestion rate
xk      <- 0.0002 # half saturation constant
hb      <- 0.17E-02   # background hazard rate (1/d)

#store parameter values at time 0 for calculating parameter changes due to stress
f<-1


km0<-km
R_m0<-R_m
f0<-1
g0<-g

#Stress function related parameters for chemical exposure---------------------------------------------------
par.kd   <- 0.482      # dominant rate constant, d-1
par.z    <- 0.112      # threshold (NEC) for survival (mg/L)
par.b   <- 1.21       # killing rate (L/mg/d) (SD only)

par.zb     <-0.0515    # no effect concentration (threshold for effect)
par.bb     <-0.908     # tolerance concentration


# Monte-Carlo repetition------------------------------------------------------------------------------------
for (m in 1:MC.no ) {
  
        food.total <- rep(c(food.daily),tmax)

        #intialize/clear state variables vectors
            ind.L        <-  c(NA)     # initial volumetric structural length (here: start at birth)
            ind.R        <-  c(NA)     # initial cumulative reproduction (here: start at birt) 
            ind.Jx       <-  c(NA)    # actual ingestion rate
            ind.e        <-  c(NA)     # scaled reserve density; e=f for constant food conditions
            ind.D      <-  c(NA)     # scaled damage
            ind.H        <-  c(NA)     # Cumulative hazard at t=0
            ind.S        <-  c(NA)     # survival probability
            ind.b_L      <-  c(NA)     # volumetric structural length of embryo
            ind.alive    <-  c(NA)     # individual still alive?
            ind.maxL     <-  c(NA)     # maximum length over time
            ind.ran      <-  c(NA)     # random number for survival
            ind.ranf     <-  c(NA)     # random number for survival
            ind.F_s      <-  c(NA)     # stress fucntion for sublethal effects
            b.L          <-  c(NA)     # structural length at egg formation
            b.R          <-  c(NA)     # calculation of brood size
            b.e          <-  c(NA)     # reserve density of egg is assumed to be constant
            ind.breeding <-  c(NA)     # does individual have embryos?


        #intial state variables of individuals in vectors
        for (i in 1:pop.ini ) {
          ind.L[i]        <-  c(Lb*sample(900:1100, 1) / 1000)    # initial length (here: start at birth)
          ind.R[i]        <-  c(0)     # initial cumulative reproduction (here: start at birt) 
          ind.Jx[i]       <-  c(0)     # actual ingestion rate
          ind.e[i]        <-  c(1)     # scaled reserve density; e=f for constant food conditions
          ind.D[i]      <-  c(0)     # scaled damage
          ind.H[i]        <-  c(0)     # Cumulative hazard at t=0
          ind.S[i]        <-  c(1)     # survival probability
          ind.alive[i]    <-  c(TRUE)  # individual still alive?
          ind.maxL[i]     <-  c(0)     # maximum length over time
          ind.ran [i]     <-  c(sample(0:1000, 1) / 1000)    # random number for survival
          ind.ranf [i]    <-  c(sample(800:1000, 1) / 1000)    # random number for survival
          ind.F_s[i]      <-  c(0)     # stress function for sublethal effects
          
          b.L[i]          <-  c(0)     # structural length at egg formation
          b.R[i]          <-  c(0)     # calculation of brood size
          b.e[i]          <-  c(1)     # reserve density of egg is assumed to be constant
          ind.breeding[i]    <-  c(0)     #  individual does not have embryos
          
          #set differnt initials for adults
          
          if (i>pop.ini-pop.adults) {
            ind.L[i]        <-  c(Lp*sample(900:1100, 1) / 1000)
            ind.H[i]        <-  c(sample(0:300, 1) / 1000)     # Cumulative hazard at t=0
            
          }
        }
        
        
        #time loop---------------------------------------------------------------------------------------------
        for (t in 1:tmax ) {
          
               #number of indivduals in vector dead or alive
               pop.size <- length(ind.L)
                 
              #---individual toxicokinetics and stress function---------------------------------------
                 for(i in 1:pop.size) {
                    #dose metric
                      ind.D[i]<-ind.D[i]+ par.kd*(Cw[t]-ind.D[i])    #scaled damage
                        #The resulting stress function for a given concentration: 
                          ind.F_s[i]  <- par.bb *max(0,ind.D[i]-par.zb)
                         
                         
                          #Parameter changes due to stress for pMoAs of fedding
                          if (pMoA == 1) {#feeding or assimilation from food
                              f  <- f0 * max(0, 1-ind.F_s[i] ) }
                     }
                     
                     
               #---population feeding--------------------------------------------  
               
                    for(i in 1:pop.size) {
                       #ind.Jx[i] <- ((sample(800:1200, 1) / 1000)* f *Jxm * ind.L[i]^ 2 * (food.total[t] / env.size)) / (xk + (food.total[t] / env.size))
                       ind.Jx[i] <- (ind.ranf[i]* f *Jxm * ind.L[i]^ 2 * (food.total[t] / env.size)) / (xk + (food.total[t] / env.size))
                       }
                 
                     # sum of feeding rates of all living individuals
                     pop.feed[t] <- sum(ind.Jx[which(ind.alive == TRUE)])
                     
                     if(pop.feed[t] > food.total[t]){
                       
                       for(i in 1:pop.size) {
                         ind.Jx[i] <- (food.total[t] / pop.feed[t]) * ind.Jx[i]
                         }
                          pop.feed[t] <- sum(ind.Jx[which(ind.alive == TRUE)])
                     }
                     
                     #add leftovers to food of next day
                     if (t+1< tmax){
                        food.total[t+1] <- food.total[t+1] + food.total[t] - pop.feed[t] }
                     
                     
                #---individual life histories--------------------------------------------
                     for(i in 1:pop.size) {      
                           if (pMoA == 2) {#somatic and maturity maintenance
                             km <- km0 * (1+ind.F_s[i] )
                             R_m <- R_m0 / (1+ind.F_s[i] )
                           }  
                           if (pMoA == 3) {#cost for structure and maturation
                             g  <- g0 * (1+ind.F_s[i] )
                             km <- km0 /(1+ind.F_s[i] )
                           } 
                           if (pMoA == 4) {#embryo hazard
                             R_m <- R_m0 * exp(-ind.F_s[i] )
                           } 
                           
                            #update scaled reserve density in iterations to minimize error
                             ind.f<-ind.Jx[i]/(Jxm*ind.L[i]^ 2)
                                for (x in 1:100 ) {ind.e[i]<-ind.e[i]+(ind.f-ind.e[i])*v/ind.L[i]/100   }
                             
                
                             # Growth
                              dL<-km*g/(3*(ind.e[i]+g)) *(ind.e[i]*v/(km*g)-ind.L[i])
                               
                               
                             #maximum length over time
                               if (ind.maxL[i] < ind.L[i]) {ind.maxL[i]<-ind.L[i] }
                              
                            
                             #Reproduction
                             if (ind.L[i]>Lp) { dR<-(R_m/((v/(km*g))^3-Lp^3)) *((v/km *ind.L[i]^2+ind.L[i]^3 )*ind.e[i]/(ind.e[i]+g)-Lp^3)}
                               else dR<-0
                            
                                 
                                 
                             #Cumulative hazard
                               ind.H[i]<-ind.H[i]+par.b*(max(0,ind.D[i]-par.z))+ hb + (s_shrink * km * max(0, ind.maxL[i] - ind.L[i])/ind.maxL[i])
                             
                             #Survival probability
                               ind.S[i]<- exp(-ind.H[i])  
                               
                               
                             #update growth and reproduction
                               ind.L[i]<-ind.L[i]+dL
                               if (dR>0) {ind.R[i]<-ind.R[i]+dR}
                            
                              #individual still alive?
                               if (ind.S[i]< ind.ran [i]) {ind.alive[i]<-FALSE  }
                               
                               
                               
                              #---calculate bood development-----------------  
                               
                               ## Reproduction
                               if ((ind.R[i] > 1) & (ind.breeding[i] == 0) & (ind.alive[i]==TRUE) )  {
                                 
                                   b.L[i] <- 0.0000001     # structural length at egg formation
                                   b.e[i] <- ind.e[i]      #reserve density of egg is assumed to be constant
                                   b.R[i]  <- trunc(ind.R[i])  # calculation of brood size
                                   ind.R[i] <-ind.R[i] -b.R[i]     # resets reproduction buffer but leaves rest of the buffer that could not be used to produce full egg
                                   
                                   ind.breeding[i] <- 1
                                 }
                                 
                                 if ((ind.breeding[i] == 1)& (ind.alive[i]==TRUE)) {
                                   b.L[i]<-b.L[i]+km*g/(3*(b.e[i]+g)) *(b.e[i]*v/(km*g)-b.L[i]) #embryo growth

                                 }
                               
                               
                               #---add newborns to population----------------------
                               if ((b.L[i]>Lb)& (ind.alive[i]==TRUE)) {
                                 for (r in 1:b.R[i] ) {  
                                     ind.breeding[i] <- 0
                                    
                                     ind.L[length(ind.L)+1]        <-  c(b.L[i])# initial length (here: start at birth)
                                     ind.R[length(ind.R)+1]        <-  c(0)     # initial cumulative reproduction (here: start at birt) 
                                     ind.Jx[length(ind.Jx)+1]       <-  c(0)     # actual ingestion rate
                                     ind.e[length(ind.e)+1]        <-  c(1)     # scaled reserve density; e=f for constant food conditions
                                     ind.D[length(ind.D)+1]      <-  c(0)     # scaled damage
                                     ind.H[length(ind.H)+1]        <-  c(0)     # Cumulative hazard at t=0
                                     ind.S[length(ind.S)+1]        <-  c(1)     # survival probability
                                     ind.alive[length(ind.alive)+1]    <-  c(TRUE)  # individual still alive?
                                     ind.maxL[length(ind.maxL)+1]     <-  c(0)     # maximum length over time
                                     ind.ran [length(ind.ran)+1]     <-  c(sample(0:1000, 1) / 1000)    # random number for survival
                                     ind.F_s[length(ind.F_s)+1]      <-  c(0)     # stress function for sublethal effects
                                     ind.ranf [length(ind.ranf)+1]    <-  c(sample(800:1000, 1) / 1000)    # random number for survival
                                     b.L[length(b.L)+1]          <-  c(0)     # structural length at egg formation
                                     b.R[length(b.R)+1]          <-  c(0)     # calculation of brood size
                                     b.e[length(b.e)+1]          <-  c(1)     # reserve density of egg is assumed to be constant
                                     ind.breeding[length(ind.breeding)+1]    <-  c(0)     # does individual not have embryos
                                 }
                                 b.R[i]<-0
                               }
                               
                               
                     }#end individual life histories
                          
                  #Store results 
                   res.pop[t,m]<-sum(ind.alive[which(ind.alive == T)])
        }  
            
  res[,1]<-c(1:tmax)
  res[,2]<-Cw
  
  mc.pop.quant <- apply(res.pop,1,quantile, na.rm = TRUE)
  
  res[,3]<-mc.pop.quant[3,] # median 
  res[,4]<-mc.pop.quant[2,] # 25%
  res[,5]<-mc.pop.quant[4,] # 75% 
  
  
}

#create plots for growth and reproduction---------------------------------------------------------
par(mfrow = c(1, 2))

barplot(res[,2], main="Exposure", ylab = "Concentration [mg/L]", xlab = "Time [d]")

#plot population sizes
plot(res [,1], res [,5], type = "l", lty=2, lwd=2, ylab = "Population size [#]", xlab = "Time [d]", main = "Population dynamics", ylim=c(0,max(res [,5])))  #
  lines(res [,1], res [,4], type = "l",lwd=2, lty=2)
   lines(res [,1], res [,3], type = "l",lwd=3)

   res.pop[70,]
        