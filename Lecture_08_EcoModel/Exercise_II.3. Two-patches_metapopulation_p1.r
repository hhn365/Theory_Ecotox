## Exercise 3. Metapopulation in two-patch system 
## Q3.1. Implement the metapopulation model in R.
## Q3.2. Why is there a decrease in the population for N2 at the beginning?
## Q3.3. What is the role of dispersal for N2 population? Test for e = 0.0 and 0.2

#############
# Load library
...


## Model function
tpdp <- new("odeModel", 
		main= function(time, init, parms) { 
			with(as.list(c(init, parms)), {
			dn1dt <- ...
			dn2dt <- ...
			list(c(dn1dt, dn2dt))
			})
		},
		parms = c(r1 = 0.2, r2 = -0.3, alpha1 = 0.01, e = 0.1, s = 0.1),
		times  = c(from = 0, to = 200, by = 1),
		init = c(n1 = 20, n2 = 20),
		solver= ...
	)	

# Model simulation
result_tpdp <- sim(...)

# Model results
...

# Model plot
# Plot two populations N1 and N2 separately
...
# Plot two populations N1 and N2 in one graph
...
