## --------------------------------------------------------------------------
library(drc)
# R package that contains the data
data(ryegrass)
head(ryegrass)


## --------------------------------------------------------------------------
ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())


## --------------------------------------------------------------------------
summary(ryegrass.LL.4)


## --------------------------------------------------------------------------
ED(ryegrass.LL.4, c(10, 50, 90))


## ----fig.width = 5, fig.height = 5-----------------------------------------
plot(ryegrass.LL.4, broken = TRUE, type = "all",
xlab = "Ferulic acid (mM)", ylab = "Root length (cm)", lwd = 2)


## ----fig.width = 5, fig.height = 5-----------------------------------------
qqnorm(residuals(ryegrass.LL.4))
abline(a = 0, b = sd(residuals(ryegrass.LL.4)))


## ----fig.width = 5, fig.height = 5-----------------------------------------
plot(fitted(ryegrass.LL.4), residuals(ryegrass.LL.4))
abline(h = 0, lty = 2)


## ----fig.width = 5, fig.height = 5-----------------------------------------
library(ggplot2)
newdata <- expand.grid(conc = exp(seq(log(0.5), log(100),
length = 100)))
# we wrap the predict function inside a function that suppresses warnings because these are non-informative in our case and would only detract in this notebook. However, the warnings would not affect functionality of the code.
suppressWarnings(pm <- predict(ryegrass.LL.4, newdata = newdata,
interval = "confidence"))
newdata$p <- pm[ , 1]
newdata$pmin <- pm[ , 2]
newdata$pmax <- pm[ , 3]
# given that the log for 0 is not defined, we set the lowest value to 0.5, but create a new concentration vector. This is good practice, never overwrite the raw data!
ryegrass$conc0 <- ryegrass$conc
ryegrass$conc0[ryegrass$conc0 == 0] <- 0.5

ggplot(ryegrass, aes(x = conc0, y = rootl)) +
geom_point() +
geom_line(data = newdata, aes(x = conc, y=p)) +
geom_ribbon(data = newdata,
aes(x = conc, y = p, ymin = pmin, ymax = pmax),
alpha = 0.2) +
coord_transform(x = "log") +
xlab("Ferulic acid (mM)") + ylab("Root length (cm)")


## --------------------------------------------------------------------------
# fit Weibull type 1 model
ryegrass.W14 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
# fit log-normal model
ryegrass.LN4 <- drm(rootl ~ conc, data = ryegrass, fct = LN.4())
# fit Brain-Cousens model
ryegrass.BC5 <- drm(rootl ~ conc, data = ryegrass, fct = BC.5())


## ----fig.width = 5, fig.height = 5-----------------------------------------
# Generating the plot
plot(ryegrass.LL.4, broken = TRUE, type="all")
# Can add the next models on top of current plot with different line types and weights
plot(ryegrass.W14, add = TRUE, col = "red", lty = 4, lwd = 1.5, broken = TRUE, type = "none")
plot(ryegrass.LN4, add = TRUE, col = "blue", lty = 2, lwd = 1.5, broken = TRUE, type = "none")
plot(ryegrass.BC5, add = TRUE, col = "darkgreen", lty = 6, lwd = 1.5, broken = TRUE, type = "none")

legend("bottomleft", legend = c("LL.4", "W1.4", "LN.4", "BC.5"),
lty = c(1, 4, 2, 6), lwd = c(1, 2, 2, 2), col = c("black", "red", "blue", "darkgreen"), bty = "n")


## --------------------------------------------------------------------------
AIC(ryegrass.LL.4)
AIC(ryegrass.W14)
AIC(ryegrass.LN4)
AIC(ryegrass.BC5)


## --------------------------------------------------------------------------
library(fitdistrplus) # to retrieve the data
data("endosulfan")
head(endosulfan)
endosulfan.art <- subset(endosulfan, group == "Arthropods" & Australian == "no")


## --------------------------------------------------------------------------
# we use the ssdtools package
library(ssdtools)  
# we fit the distribution
ssd_ln <- ssd_fit_dists(data = endosulfan.art, left = "ATV", dists = "lnorm")
# and create the plot
ssd_plot_cdf(ssd_ln)


## ----fig.width = 5, fig.height = 5-----------------------------------------
# First we extract and sort the concentrations
Conc <- sort(endosulfan.art$ATV)  
# Extract estimated parameters
params <- ssd_ln$lnorm$pars
# Assign mean and standard deviation to new object
# But need to convert to object type numeric
mlog <- as.numeric(params["meanlog"])
sdlog <- as.numeric(params["log_sdlog"])
# Now we determine the theoretical log-normal distribution
# Note that you would need to adjust this, if you fitted a different distribution
q_dist <- qnorm(ppoints(length(Conc)), mean = mlog, sd = sdlog)

# Generate Q-Q plot
qqplot(q_dist, log10(Conc), xlab = "Theoretical quantiles", ylab = "Observed quantiles",
  main = "QQ-Plot (lognormal fit)")
# We provide the ordinates for probability plotting with ppoints, the center of the distribution and the variability of the distribution. Note that once we know the center and the variability along with the sample size, we can specify the exact values of a theoretical distribution.
qqline(log10(Conc), distribution = function(p) qnorm(p, mean = mlog, sd = sdlog), col = "red")


## --------------------------------------------------------------------------
ssd_hc(ssd_ln, proportion = 0.05)
# we can also generate bootstrapped confidence intervals for this proportion
ssd_hc(ssd_ln, proportion = 0.05, ci = TRUE)


## --------------------------------------------------------------------------
ssd_fitdists <- ssd_fit_dists(data = endosulfan.art, left = "ATV", dists = c("lnorm", "burrIII3", "weibull", "llogis"))
ssd_plot_cdf(ssd_fitdists)


## --------------------------------------------------------------------------
ssd_gof(ssd_fitdists, wt = TRUE)


## --------------------------------------------------------------------------
hc_table <- ssd_hc(ssd_fitdists, proportion = 0.05, average = FALSE, ci = TRUE)
print(hc_table)
# let us also plot the estimates with the confidence intervals.
library(ggplot2)
ggplot(hc_table, aes(x = dist, y = est)) +
  geom_point(size = 3) +  # Plot estimates
  geom_errorbar(aes(ymin = lcl, ymax = ucl), width = 0.2) +  # Confidence intervals
  labs(
    x = "Distribution",
    y = "HC5",
    title = "Estimates with 95% Confidence Intervals"
  ) +
  theme_minimal()




## --------------------------------------------------------------------------
library(drc)
# Load data set
data(S.alba)
# Brief look at the data
head(S.alba)
# 
salba_ll4 <- drm(DryMatter ~ Dose, curveid = Herbicide, data = S.alba, fct = LL.4())
# Note: If prior information is available regarding the similarity of model parameters such as the upper or lower limit, it can be incorporated using the argument pmodel.

# Display the model results 
summary(salba_ll4)



## ----fig.width = 5, fig.height = 5-----------------------------------------
plot(salba_ll4, broken = TRUE, xlab = "Dose (g a.i./ha)", ylab = "Dry matter (g/pot)")


## --------------------------------------------------------------------------
EDcomp(salba_ll4, percVec = c(50, 50), interval = "delta")


## --------------------------------------------------------------------------
compParm(salba_ll4, strVal = "e", operator = "-")
# the operator determines the comparison type. Here, "-" tests the difference. Use "/" to compare parameter ratios.
compParm(salba_ll4, strVal = "b", operator = "-")


## ----eval=FALSE------------------------------------------------------------
# # This code will be shown but not executed
# # You need to supply your object before the code works
# library(multcomp)
# summary(glht(<object>[["EDmultcomp"]], linfct = contrMat(1:4, "Tukey")))
# # This performs pairwise comparisons using Tukey's HSD test, which is commonly used as post-hoc test after an ANOVA.

