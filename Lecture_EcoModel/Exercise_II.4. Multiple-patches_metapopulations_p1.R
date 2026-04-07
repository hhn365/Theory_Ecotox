# Exercise 4. Metapopulation expansion over a generic landscape
# Reference: Mestre et al. 2016 Environ. Model. Softw. 81: 40-44
################################################################

# Load library
library(MetaLandSim)

######################
## Generate landscape 
set.seed(123)
rland1 <- rland.graph(mapsize = ,  
                   dist_m = , 
                   areaM = , 
                   areaSD = , 
                   ..., 
                   ...,    
                   plotG = TRUE)

# Landscape summary
class(...)
summary_landscape(...)

# Node characteristics
nodes <- rland1$nodes.characteristics
nodes
text(x = nodes[,"x"],y = nodes[,"y"], pos = 2, 
     labels = as.character(nodes[,"ID"]))

# Edge characteristics
e1 <- edge.graph(...) 
e1

######################
## Simulate metapopulation occupancy of landscape
set.seed(...)
spec1 <- species.graph(rl = ..., 
                       ..., 
                       ..., 
                       nsew = "none", 
                       plotG = TRUE) 
# label the node points
text(x = nodes[,"x"],y = nodes[,"y"], pos = 2, 
     labels = as.character(nodes[,"ID"]))

class(...)
summary_metapopulation(...)
spec1$nodes.characteristics

######################
## Simulate metapopulation range expansion
# Species span
span1 <- span.graph(rl=..., 
                    span=..., 
                    par1="none")

# Model parameters
param1 <- create.parameter.df(alpha=..., ..., ..., ...)

# Model simulation
sim1 <- simulate_graph(rl=...,
                      rlist=...,
                      simulate.start=FALSE,
                      nsew="none",
                      succ = "none",
                      param_df=param1,
                      kern="op1",
                      conn="op1",
                      colnz="op1",
                      ext="op1",
                      b=1)

# Model plot 
par(mfrow=c(1,3))

plotL.graph(rl=spec1, rlist=sim1, 
nr=1, species=TRUE, links=FALSE)

plotL.graph(rl=spec1, rlist=sim1, 
nr = 20, ..., ...)

plotL.graph(rl=spec1, rlist=sim1, 
nr = 100, ..., ...)

par(mfrow = c(1, 1)) 
