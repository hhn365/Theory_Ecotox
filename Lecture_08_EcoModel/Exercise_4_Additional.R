# References: Codes are adapted from: 
# 1. Stevens 2009 - A Primer of Ecology with R
# 2. Introduction to theoretical ecology (https://pojuke.github.io/TheoreticalEcologyPJK/week-8---lotka-volterra-competition-model---population-dynamics.html)

########################
# Load library
library(ggplot2)
library(tidyverse)
library(deSolve)
library(grid)  

################## Stability analysis - Case Unstable coexistence##########
## 1. Assign values of model parameters (Jacobian analysis)
r1 <- 1; r2 <- 0.1
K1 <- 100; K2 <- 50
a11 <- 1 / K1       # intraspecific competition for species 1
a22 <- 1 / K2       # intraspecific competition for species 2
a12 <- 0.2 / K1     # interspecific competition effect of N2 on N1
a21 <- 0.5 / K2     # interspecific competition effect of N1 on N2

########################
## 2. Model expression for dN/dt
dN1dt <- expression(r1*N1 - r1*a11*N1^2 - r1*a12*N1*N2)
dN2dt <- expression(r2*N2 - r2*a22*N2^2 - r2*a21*N1*N2)

# Partial derivatives (Jacobian components)
(ddN1dN1 <- D(dN1dt, "N1"))
(ddN1dN2 <- D(dN1dt, "N2"))
(ddN2dN1 <- D(dN2dt, "N1"))
(ddN2dN2 <- D(dN2dt, "N2"))

# Jacobian matrix 
J <- expression(matrix(
  c(eval(ddN1dN1), eval(ddN1dN2),
    eval(ddN2dN1), eval(ddN2dN2)),
  nrow = 2, byrow = TRUE))

########################
## 3. Calculate model equilibrium (at N1*, N2*)
N1Star <- expression((a22 - a12) / (a22*a11 - a12*a21))
N2Star <- expression((a11 - a21) / (a22*a11 - a12*a21))

N1 <- eval(N1Star)
N2 <- eval(N2Star) # N2* is zero -> no coexistence!

########################
## 4. Evaluate Jacobian matrix at equilibrium
(J1 <- eval(J))
eigStable <- eigen(J1)
eigStable[["values"]] 

########################
## 5. Plotting function - Phase plane
phase_plane <- function(r1, r2, a11, a21, a22, a12, title, t){
  # Vectors
  LV <- function(times, state, parms) {
    with(as.list(c(state, parms)), {
      dN1_dt = N1 * (r1 - a11*N1 - a12*N2)
      dN2_dt = N2 * (r2 - a22*N2 - a21*N1)
      return(list(c(dN1_dt, dN2_dt)))
    })
  }
  # Parameters
  times <- c(0, t)
  parms <- c(r1 = r1, r2 = r2, a11 = a11, a21 = a21, a22 = a22, a12 = a12)
  # Plot range with a buffer expansion for better visualization
  x_inter<- max(c(r1/a11, r2/a21)) *1.2
  y_inter <- max(c(r2/a22, r1/a12)) *1.2
  # Arrows
  vector_grid <- expand.grid(seq(5, x_inter, length.out = 10),
                             seq(5, y_inter, length.out = 10))
  vector_data <- vector_grid %>%
    pmap(., function(Var1, Var2){
      state <- c(N1 = Var1, N2 = Var2)
      pop_size <- ode(func = LV, times = times, y = state, parms = parms)
      pop_size[2, 2:3]
    }) %>%
    bind_rows() %>%
    rename(xend = N1, yend = N2) %>%
    bind_cols(vector_grid) %>%
    rename(x = Var1, y = Var2)
    # Phase plane
    ggplot() +
      geom_abline(slope = -a11/a12, intercept = r1/a12, 
                  color = "orange", size = 1.5) +
      geom_abline(slope = -a21/a22, intercept = r2/a22, 
                  color = "blue", size = 1.5) +
      geom_segment(data = vector_data,
                   aes(x = x, y = y, xend = xend, yend = yend),
                   arrow = arrow(length = unit(0.1, "cm"))) +    
      scale_x_continuous(name = "N1", limits = c(0, x_inter), expand = c(0, 0)) +
      scale_y_continuous(name = "N2", limits = c(0, y_inter), expand = c(0, 0)) +
      theme_bw(base_size = 13) +
      theme(panel.grid = element_blank(),
            plot.title = element_text(hjust = 0.5),
            aspect.ratio = 1) +
      labs(title = title)
  }

########################
## 6. Phase plane plot
phase_plane(r1 = r1, r2 = r2, a11 = a11, a21 = a21,
            a22 = a22, a12 = a12, t = 0.2,
            title = "Phase Plane - Competitive exclusion")

################## Stability analysis - Case Stable coexistence ##########
# Assign new values
a11 <- a22 <- 0.01
a12 <- a21 <- 0.001 
r1 <- r2 <- 1

# Values in equilibrium for N
N1Star <- expression((a22-a12)/(a22*a11 - a12*a21))
N2Star <- expression((a11-a21)/(a22*a11 - a12*a21))
N1 <- eval(N1Star)
N2 <- eval(N2Star)
# N1 and N2 > zero - coexistence!

# Evaluate Jacobian matrix
(J2 <- eval(J))
eigStable <- eigen(J2); eigStable[["values"]] 

phase_plane(r1 = r1, r2 = r2, a11 = a11, a21 = a21,
            a22 = a22, a12 = a12, t = 0.2,
            title = "Phase Plane - Stable coexistence")
# Comment: All real parts < 0 -> globally stable point and no imagery parts (i.e., no oscillation)

################## Stability analysis - Case N2 win ##########
# Assign new values
a11 <- a22 <- 0.01
a12 <- 0.01 
a21 <- 0.001
r1 <- r2 <- 1

# Values in equilibrium for N
N1Star <- expression((a22-a12)/(a22*a11 - a12*a21))
N2Star <- expression((a11-a21)/(a22*a11 - a12*a21))
N1 <- eval(N1Star)
N2 <- eval(N2Star)
# N1 = 0 - no coexistence!

# Evaluate Jacobian matrix
(J2 <- eval(J))
eigStable <- eigen(J2); eigStable[["values"]] 

phase_plane(r1 = r1, r2 = r2, a11 = a11, a21 = a21,
            a22 = a22, a12 = a12, t = 0.2,
            title = "Phase Plane - Competitive exclusion")
