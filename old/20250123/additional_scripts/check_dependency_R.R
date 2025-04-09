#! /bin/env Rscript


#####################################################################
# List of packages to check
packages_to_check <- c("getopt", "parallel", "matrixStats", "seqinr", "phytools", "phangorn", "simclock1")

# Get a list of all installed packages
installed_packages <- installed.packages()[, "Package"]

# Check if each package is installed
are_installed <- sapply(packages_to_check, function(pkg) pkg %in% installed_packages)

# Print the results
for (pkg in packages_to_check) {
  if (are_installed[pkg]) {
    #cat(pkg, "is installed.\n")
  } else {
    cat("R package ", pkg, " is NOT installed. Try install.packages('", pkg, "\')", " or devtools::install_github('dosreislab/simclock')", "\n", sep = "")
  }
}

