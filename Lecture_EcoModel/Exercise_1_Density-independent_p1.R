# Question 1- Exponential population growth model

# Load library
library(deSolve)

# Model function - Exponential growth dN/dt = r*N
dens_indep <- function(times, init, parms){
  with(as.list(c(init,parms)), {
    dN_dt = r
    return(list(c(dN_dt)))
  })
}

# Model parameters
times   <- seq(from = 0, to = 120, by = 1)
initial <- c(N = 10)
parms   <- c(r = 0.03)

# Model simulation run using ODE
result_dens_indep <- ode(func = dens_indep, times = times, y = init)
head(result_dens_indep)

# Model plot
plot(N ~ time, data = result_dens_indep, pch = 1, cex = 0.8, main = "Density-independent growth model")

# Optional- Analytical solution: N(t) = N_0 * e^(rt)
curve(init[1]*exp(parms[1]*x), col = "red", lwd = 2, add = T)

## Questions:
# Question_1. Run the model in r - resolve errors in model function, model parameters, and model simulation 
# Question_2. Test if changing population density N and growth rate r affect the result pattern


