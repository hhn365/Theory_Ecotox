
## Important notice - read before you run the code 
## You need to manually install JAGS on your computer: https://sourceforge.net/projects/mcmc-jags/files/
## This means you need to download and install locally
## If you do not have JAGS installed, R packages will not work
## You also should have installed JAVA, if this is not the case, you can install it here
## again manually:
## https://www.java.com/de/download/

pkgs <- readLines(file("https://raw.githubusercontent.com/hhn365/Theory_Ecotox/refs/heads/main/Lecture_Introduction/installed_pkgs.txt", "r"))
str(pkgs)
install.packages(pkgs)
BiocManager::install("DESeq2")
install.packages("remotes")
remotes::install_github('andschar/standartox')
remotes::install_github('andschar/chemlook')
remotes::install_github("alexology/biomonitoR", ref = "main", build_vignettes = TRUE)