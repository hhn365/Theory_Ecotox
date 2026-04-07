## Exercise  1. Levins metapopulation model
## Q1.Implement the metapopulation model in R. 
##    Define the equilibrium.

#############
# Load library
...

#############
## Q1. Levins metapopulation dp_dt = ci*p*(1-p) - e*p
## Model function 
levins <- function(t, init, parms) {
       p <- init[1]            # init is the initial state vector
       with(as.list(parms), {
       dp <- ...
       return(list(...))
       }) 
}

# Model parameters
times <- ...
...
...

# Model simulation using ODE
out.L <- ode(y=init, func=..., ..., ...)
head(...)   # read results

# Model plot
plot(x = out.L[,1], y = ..., 
     type = "l",col = "blue", lwd = 2,                    
     ylim = c(0, 1), 
     xlab = "time t", ylab = ...,
     main = "Levins model")

text(35, 0.60, labels = bquote(frac(dp, dt) == cp*(1 - p) - ep))

# Model equilibrium: dN/dt = 0
...
