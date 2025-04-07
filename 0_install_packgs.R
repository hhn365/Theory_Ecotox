
## Important notice - read before you run the code 
## You need to manually install JAGS on your computer: https://sourceforge.net/projects/mcmc-jags/files/
## This means you need to download and install locally
## If you do not have JAGS installed, R packages will not work

pkgs <- readLines(file("https://raw.githubusercontent.com/hhn365/Theory_Ecotox/refs/heads/main/Lecture_01_Introduction/installed_pkgs.txt", "r"))
str(pkgs)
install.packages(pkgs)
BiocManager::install("DESeq2")
