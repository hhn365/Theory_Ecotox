# Load library
library(simecol)

## Question 1 - Implement the species competition model
lvcomp2 <- new("odeModel", 
		main= function(time, init, parms) { 
			with(as.list(c(init, parms)), {
			dn1dt <- r1 * n1 * (1 - ((n1 + a12 * n2)/K1)) 
			dn2dt <- r2 * n2 * (1 - ((n2 + a21 * n1)/K2)) 
			list(c(dn1dt, dn2dt))
			})
		},
		parms = c(r1 = 1, r2 = 0.1, a12 = 0.2, a21 = 0.5,K1 = 100, K2 = 50),
		times = c(from = 0, to = 100, by = 1),
		init = c(n1 = 20, n2 = 20),
		solver= "lsoda"
	)	

# Run simulation 
res_1 <- sim(lvcomp2)

# Model results
print(res_1, all=TRUE)
head(res_1@out) 

## Model plots
plot(res_1, las = 1, lwd = 2, col = "blue")
par(mfrow = c(1, 1)) # Reset to single plot layout

plot(res_1@out[,"time"], res_1@out[,"n1"], 
     type = "l", col = "blue", lty=2, lwd = 2,
     xlab = "Time", ylab = "N1 and N2",
     ylim = range(res_1@out[,c("n1", "n2")]),
     main = "Competition Model")
lines(res_1@out[,"time"], res_1@out[,"n2"], col = "orange", lwd = 2)
legend(80,70, legend = c("N1", "N2"),
       col = c("blue", "orange"), lwd = 2)

####################
## Question 2 - Find out for which values of K and a the species can coexist  					
#  and under which conditions only one survives			


