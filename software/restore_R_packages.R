# restore_packages.R
#
# installs each package from the stored list of packages

load("/tmp/installed_packages.rda")
sink("Rinstallation.output")

if (!requireNamespace("BiocManager", quietly = TRUE)) #it will not be installed
    install.packages("BiocManager")
BiocManager::install() #optional to install core packages

# this will install R CRAN packages
install.packages(installedpackages,repos = "http://cran.us.r-project.org")
update.packages()

# this will install Bioconductor packages
for (i in 1:length(installedpackages)) BiocManager::install(installedpackages[i])

sink()
