## Exercise  1. Levins metapopulation model
## Q1.Implement the metapopulation model in R. 
##    Define the equilibrium.
## Q2. Compare Levins model with population density-dependent growth model.

#############
# Load library
library(deSolve)

#############
## Q1. Levins metapopulation dp_dt = ci*p*(1-p) - e*p
## Model function 
levins <- function(t, init, parms) {
       p <- init[1]            # init is the initial state vector
       with(as.list(parms), {
       dp <- ci * p * (1-p) - e * p
       return(list(dp))
       }) 
}

# Model parameters
times <- seq(from = 1, to = 100, by = 1)
parms <- c(ci=0.15, e=0.05)
init <- 0.01

# Model simulation using ODE
out.L <- ode(y=init, times=times, func=levins, parms=parms)
head(out.L)

# Model plot
plot(x = out.L[,1], y = out.L[,2], 
     type = "l",col = "blue", lwd = 2,                    
     ylim = c(0, 1), 
     xlab = "time t", ylab = "proportion of patches p",
     main = "Levins model")
text(35, 0.60, labels = bquote(frac(dp, dt) == cp*(1 - p) - ep))

# Model equilibrium: dN/dt = 0
p_star_L <- 1- parms[2]/parms[1]

## Q2. Compare Levins model with population density-dependent growth model
# Solution: in the lecture presentation (slide 11)
