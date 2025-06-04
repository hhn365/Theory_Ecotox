# Load library
library(simecol) # Simulation of Ecological Dynamic Systems
# install.packages("simecol")

## Q1. Run the density-dependent population growth model 
# Model using differential equation ODE 
dens_dep <- new("odeModel",                # model class  (odeModel, gridModel, rwalkModel, indbasedModel, etc.)  
		main = function(times, init, parms) { 
			with(as.list(c(init, parms)), {
			dN_dt <- r * N * (1-alpha*N)
			list(c(dN_dt))
			})
		},
		parms  = c(r = 0.3, alpha = 0.005),
		times  = c(from = 0, to = 120, by = 1),
		init   = c(N=2),
		solver= "lsoda"  # numerical algorithm
			)	
# Run simulation
result_ddep <- sim(dens_dep)

# Model results
print(result_ddep, all=TRUE)  # print all results
head(result_ddep@out)         # show only results from output matrix 'out'

# Model plot
plot(result_ddep@out[ , 1], result_ddep@out[ , 2], type = "l", xlab = "Time", ylab = "Count (N)", 
     las = 1, lwd = 2, col = "blue")

########################
## Question 2 - Manipulate r and α and see how the model reacts

## Question 3 - For which values of N is the model in equilibrium 
# Note: Model in equilibrum when dN/dt = 0

## Question 4 - Which equilibrium is stable?

######################## Example solutions #############
## Solution to question 2
# Test different parameter sets using a loop
r_val <- c(0.2, 0.3, 0.4)        # r = 0.3 
alpha_val <- c(0.005, 0.01)      # alpha = 0.005

# Plot setup
plot(NULL, xlim = c(0, 100), ylim = c(0, 200), xlab = "Time", ylab = "Count (N)",
     las = 1, main = "Population Growth under Different r and alpha")

# Color palette for plotting
colors <- rainbow(length(r_val) * length(alpha_val))
legend_labels <- c()

# Loops
i <- 1
for (r in r_val) {
  for (alpha in alpha_val) {
    # Copy the model and test new parameters
    dens_dep <- dens_dep
    parms(dens_dep) <- c(r = r, alpha = alpha)
    # Run models with new parameters
    result <- sim(dens_dep)
    out <- as.data.frame(result@out)
    # Plot the result
    lines(out$time, out$N, col = colors[i], lwd = 2)
    legend_labels[i] <- paste0("α = ", alpha, ", r = ", r)
    i <- i + 1
    }
}
legend("bottomright", legend = legend_labels, col = colors, lwd = 1, cex = 0.8, bty = "n")

## Solution to question 3
# Model is in equilibrium in two cases: N = 0 or N = 1/alpha (non-zero equilibrium)
alpha = 0.005
N_equilibrium <- c(0, 1/alpha)

## Solution to question 4 - Stability at equilibrium 
# Plot the growth rate against the population size 
# Plot setting
par(mar = c(4, 4, 2, 1), cex = 1.2, las = 1)

# Population sizes and expression
r <- 0.3
alpha <- 0.005
N <- 0:270   # N* = 1/alpha = 1/0.005 = 200
pop.growth.rate <- expression(r * N * (1 - alpha*N))

# Plot curve
plot(N, eval(pop.growth.rate), type = "l", , lwd = 2, col = "blue",
     ylab = "dN/dt", xlab = "N",
     main = "Population growth rate")

# Equilibrium line
abline(h = 0)   

# Legend
legend('topright', "r = 0.3", lty = 1, lwd = 2, cex = 0.8)

# Points
N <- c(0, 50, 150, 200, 250)
points(N, eval(pop.growth.rate), cex = 1.5)
text(N, eval(pop.growth.rate), letters[1:5],adj = c(.5, 2))

# Arrows
arrows(140, 4, 170, 4, length = .1, lwd = 3)
arrows(250, -4, 220, -4, length = .1, lwd = 3)

# dev.print(pdf, "Stability1.pdf")