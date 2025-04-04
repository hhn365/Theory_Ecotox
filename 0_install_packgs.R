pkgs <- readLines(file("https://raw.githubusercontent.com/hhn365/Theory_Ecotox/refs/heads/main/Lecture_01_Introduction/installed_pkgs.txt", "r"))
str(pkgs)
install.packages(pkgs)
BiocManager::install("DESeq2")

## You need to manually install JAGS on your computer: https://sourceforge.net/projects/mcmc-jags/files/
## This means you need to download and install locally 
