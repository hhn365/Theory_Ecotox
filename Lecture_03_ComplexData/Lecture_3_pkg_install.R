# Install packages

# Option 1 - From Lecture 1
pkgs <- readLines(file("https://raw.githubusercontent.com/hhn365/Theory_Ecotox/refs/heads/main/Lecture_03_ComplexData/installed_pkgs.txt", "r"))
str(pkgs)
install.packages(pkgs)

# Option 2
pkgs <- c("ggplot2", "dplyr", "lattice", "MASS", "lme4", "nlme")
install.packages(pkgs)
