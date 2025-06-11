# Exercise 4. Metapopulation expansion over a generic landscape
# Reference: Mestre et al. 2016 Environ. Model. Softw. 81: 40-44
################################################################

# Load library
library(MetaLandSim)

######################
## Generate landscape 
set.seed(123)
rland1 <- rland.graph(mapsize = 1000,  
                   dist_m = 100, 
                   areaM = 0.2, 
                   areaSD = 0.1, 
                   Npatch = 20, 
                   disp = 150,    
                   plotG = TRUE)

# Landscape summary
class(rland1)
summary_landscape(rland1)

# Node characteristics
nodes <- rland1$nodes.characteristics
nodes
text(x = nodes[,"x"],y = nodes[,"y"], pos = 2, labels = as.character(nodes[,"ID"]))

# Edge characteristics
e1 <- edge.graph(rland1) 
e1

######################
## Simulate metapopulation occupancy of landscape
set.seed(124)
spec1 <- species.graph(rl = rland1, 
                       method = "percentage", 
                       parm = 20, 
                       nsew = "none", 
                       plotG = TRUE) 
# label the node points
text(x = nodes[,"x"],y = nodes[,"y"], pos = 2, labels = as.character(nodes[,"ID"]))

class(spec1)
summary_metapopulation(spec1)
spec1$nodes.characteristics

######################
## Simulate metapopulation range expansion
# Species span
span1 <- span.graph(rl=rland1, span=100, par1="none")

# Model parameters
param1 <- create.parameter.df(alpha=0.003,  # inverse dispersal distance
                              x=0.5,        # Scaling exponent of extinction-area (larger x will decrease extinction risk as patch size increase)
                              y=2,          # Scaling component for colonization probability          
                              e=0.01)       # Global extinction probability of all patches (increase e will increase extinction in all patches)

# Model simulation
sim1 <- simulate_graph(rl=spec1,
                      rlist=span1,
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
nr=20, species=TRUE, links=FALSE)

plotL.graph(rl=spec1, rlist=sim1, 
nr=100, species=TRUE, links=FALSE)

par(mfrow = c(1, 1)) 
