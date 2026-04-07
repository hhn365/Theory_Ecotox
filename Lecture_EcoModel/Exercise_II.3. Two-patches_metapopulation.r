## Exercise 3. Metapopulation in two-patch system 
## Q3.1. Implement the metapopulation model in R.
## Q3.2. Why is there a decrease in the population for N2 at the beginning?
## Q3.3. What is the role of dispersal for N2 population? Test for e = 0.0 and 0.2

#############
# Load library
library(simecol)

## Model function
tpdp <- new("odeModel", 
		main= function(time, init, parms) { 
			with(as.list(c(init, parms)), {
			dn1dt <- r1 * n1 * (1-alpha1*n1) - e * n1 * (alpha1*n1)^s  + e *n2
			dn2dt <- r2 * n2                 + e * n1 * (alpha1*n1)^s  - e * n2  
			list(c(dn1dt, dn2dt))
			})
		},
		parms = c(r1 = 0.2, r2 = -0.3, alpha1 = 0.01, e = 0.1, s = 0.1),
		times  = c(from = 0, to = 200, by = 1),
		init = c(n1 = 20, n2 = 20),
		solver= "lsoda"
	)	

# Model simulation
result_tpdp <- sim(tpdp)

# Model results
print(result_tpdp, all=TRUE)
head(result_tpdp@out)

# Model plot
plot(result_tpdp, las = 1, lwd = 2, col = "blue")

# Plot both populations in same graph
par(mfrow = c(1, 1)) 

plot(result_tpdp@out[,"time"], result_tpdp@out[,"n1"], 
     type = "l", col = "blue", lty=2, lwd = 2,
     xlab = "Time", ylab = "N1 and N2",
     ylim = range(result_tpdp@out[,c("n1", "n2")]),
     main = "Two patch metapopulation model")
lines(result_tpdp@out[,"time"], result_tpdp@out[,"n2"], col = "orange", lwd = 2)
legend(150,50, legend = c("N1", "N2"), col = c("blue", "orange"), lwd = 2)

##############################
# Solution - Q3.2. 
# Reduction because of negative growth rate r2 = -0.3 and N2 would go extinct
# Constant migration leads to stable population

##############################
# Solution - Q3.3. 
# e = 0 equal no emigration of population from N1 to N2 => N2 goes extinct
# e = 0.2 > 0.1, higher population emigration from N1 to N2 => N1 decreases and N2 increases

# Model function
tpdp2 <- new("odeModel", 
		main= function(time, init, parms) { 
			with(as.list(c(init, parms)), {
			dn1dt <- r1 * n1 * (1-alpha1*n1) - e * n1 * (alpha1*n1)^s  + e *n2
			dn2dt <- r2 * n2                 + e * n1 * (alpha1*n1)^s  - e * n2  
			list(c(dn1dt, dn2dt))
			})
		},
		parms = c(r1 = 0.2, r2 = -0.3, alpha1 = 0.01, e = 0.2, s = 0.1),
		times  = c(from = 0, to = 200, by = 1),
		init = c(n1 = 20, n2 = 20),
		solver= "lsoda"
	)	

# Model simulation
result_tpdp2 <- sim(tpdp2)

# Model results
print(result_tpdp2, all=TRUE)
head(result_tpdp2@out)

# Model plot
plot(result_tpdp2, las = 1, lwd = 2, col = "blue")

# Plot both populations in same graph
par(mfrow = c(1, 1)) 

# Plot compare N2 in model with s = 0.1 and S = 0.5
plot(result_tpdp@out[,"time"], result_tpdp@out[,"n1"], 
     type = "l", col = "blue", lty=2, lwd = 2,
     xlab = "Time", ylab = "N1 and N2",
     #ylim = range(result_tpdp2@out[,c("n1", "n2")]),
     ylim = c(0, 220),
     main = "Two-patch model: e = 0.0")
lines(result_tpdp@out[,"time"], result_tpdp@out[,"n2"], col = "orange", lwd = 2)
lines(result_tpdp2@out[,"time"], result_tpdp2@out[,"n1"], col = "lightblue", lty=2, lwd = 2)
lines(result_tpdp2@out[,"time"], result_tpdp2@out[,"n2"], col = "red", lwd = 2)
legend(40,220, legend = c("N1", "N1 - no dispersal", "N2", "N2 - no dispersal"), col = c("blue","lightblue", "orange", "red"), lwd = 2)
