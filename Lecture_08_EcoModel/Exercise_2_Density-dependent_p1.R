# Load library
library(simecol) # Simulation of Ecological Dynamic Systems
                 # install.packages("simecol")

## Question_1. Run the density-dependent population growth model 
# Model using differential equation ODE 
dens_dep <- new("odeModel",                # model class  (odeModel, gridModel, rwalkModel, indbasedModel, etc.)  
		main = function(times, init, parms) { 
			with(as.list(c(init, parms)), {
			dN_dt <- r * N * (1-alpha*N)
			list(c(dN_dt))
			})
		},
		parms = c(),     # provide input data for the variable
		times = c(),     # provide input data for the variable 
		init = c(),      # provide input data for the variable
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
## Question_2 - Manipulate r and α and see how the model reacts

## Question_3 - For which values of N is the model in equilibrium 

## Question_4 - Which equilibrium is stable?
