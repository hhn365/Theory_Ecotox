## Exercise 2. Implement three model types in R.
## Q2.1. Implement three metapopulation models: Levins, Propagule rain, Core-satellite.
## Q2.2. Plot them jointly in one plot for an arbitrary set of parameters.
## Q2.3. Determine equilibrium conditions of each model (dp/dt =0).

#############
# Load library
library(deSolve)

#############
# Model parameters
times <- seq(from = 1, to = 100, by = 1)
init <- 0.01
e <- 0.05
ci <- ce <- c <- 0.15
parms <- c(c=0.15, e=0.05)

#############
# Model functions 
# Model 1 - Levins metapopulation model 
levins <- function(t, init, parms) {
       p <- init[1]            
       with(as.list(parms), {
       dp <- ci * p * (1-p) - e * p
       return(list(dp))
       }) 
}

# Model 2 - Propagule rain metapopulation model 
gotelli <- function(t, init, parms) {
       p <- init[1]            
       with(as.list(parms), {
       dp <- ce * (1-p) - e * p
       return(list(dp))
       }) 
}

## Model 3 - Core-satellite metapopulation model
hanski <- function(t, init, parms) {
       p <- init[1]            
       with(as.list(parms), {
       dp <- ci * p * (1-p) - e * p * (1-p)
       return(list(dp))
       }) 
}

#############
# Model simulation using ODE
out.L <- data.frame(ode(y=init, times=times, func=levins, parms=parms))
out.G <- data.frame(ode(y=init, times=times, func=gotelli, parms=parms))
out.H <- data.frame(ode(y=init, times=times, func=hanski, parms=parms))

head(out.L) # read results
head(out.G)
head(out.H)

#####################
# Model plots - Joint graph
plot(out.L[,1], out.L[,2], 
     type = "l", col = "blue", lwd = 2,
     ylim = c(0, 1), 
     xlab = "time t", ylab = "proportion of patches p",
     main = "Metapopulation Models")
text(47, 0.4, labels = expression("Levins: " ~ frac(dp, dt) == c[i]*p*(1 - p) - e*p), col = "blue", pos = 4)

lines(out.G[,1], out.G[,2], col = "darkgreen", lwd = 2)
text(15, 0.85, labels = expression("Propagule rain: " ~ frac(dp, dt) == c[e]*(1 - p) - e*p), col = "darkgreen", pos = 4)

lines(out.H[,1], out.H[,2], col = "orange", lwd = 2)
text(24, 0.6, labels = expression("Core-satellite: " ~ frac(dp, dt) == c[i]*p*(1 - p) - e*p*(1 - p)), col = "orange", pos = 4)

#####################
## Determine equilibrium conditions of each model (dp/dt =0)
# Levins: dp/dt = ci * p * (1 - p) - e * p
(p_star_L <- 1- e/ci)
# Propagule rain: dp/dt = ce * (1 - p) - e * p
(p_star_P <- ce/(ce + e))
# Core-satellite: dp/dt = ci * p * (1 - p) - e * p * (1 - p)
(p_star_H1 <- 1)
(p_star_H2 <- 0)

