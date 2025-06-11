## Exercise 2. Implement three model types in R.
## Q2.1. Implement three metapopulation models: Levins, Propagule rain, Core-satellite.
## Q2.2. Plot them jointly in one plot for an arbitrary set of parameters.
## Q2.3. Determine equilibrium conditions of each model (dp/dt =0).

#############
# Load library
...

#############
# Model parameters
times <- seq(...)
init <- ...
e <- ...
ci <- ce <- c <- ...
parms <- ...

#############
# Model functions 
# Model 1 - Levins metapopulation model 
levins <- function(t, init, parms) {
       p <- init[1]            
       with(as.list(parms), {
       dp <- ...
       return(list(dp))
       }) 
}

# Model 2 - Propagule rain metapopulation model 
gotelli <- function(t, init, parms) {
       p <- init[1]           
       with(as.list(parms), {
       dp <- ...
       return(list(...))
       }) 
}

## Model 3 - Core-satellite metapopulation model
hanski <- function(t, init, parms) {
       ...            
       ...
       ...
       ...
       }) 
}

#############
# Model simulation using ODE
out.L <- data.frame(ode(y=init, times=times, func=..., parms=parms))
out.G <- data.frame(ode(..., ..., ..., ...))
out.H <- data.frame(ode(..., ..., ..., ...))

head(...) # read results
head(...)
head(...)

#############
# Model plots - Joint graph
# Levins
plot(out.L[,1], ..., 
     type = "l", col = "blue", lwd = 2,
     ylim = c(0, 1), 
     xlab = "time t", ylab = ...,
     main = "Metapopulation Models")
text(47, 0.4, labels = expression("Levins: " ~ frac(dp, dt) == c[i]*p*(1 - p) - e*p), col = "blue", pos = 4)

# Propagule rain
lines(..., ..., col = "darkgreen", lwd = 2)
text(15, 0.85, labels = expression("Propagule rain: " ~ frac(dp, dt) == ...), col = "darkgreen", pos = 4)

# Core-satellite
lines(..., ..., col = "orange", lwd = 2)
text(24, 0.6, labels = expression("Core-satellite: " ~ frac(dp, dt) == ...), col = "orange", pos = 4)

################
# Determine equilibrium conditions of each model (dp/dt =0)
# Levins: 
(p_star_L <- ...)
# Propagule rain: 
(p_star_P <- ...)
# Core-satellite: 
(p_star_H <- ...)

