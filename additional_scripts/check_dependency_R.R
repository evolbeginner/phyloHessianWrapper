#! /bin/env Rscript


#####################################################################
options(repos = c(CRAN = "https://cloud.r-project.org"))

args <- commandArgs(trailingOnly = TRUE)


#####################################################################
# List of packages to check
packages_to_check <- c("getopt", "parallel", "seqinr", "phytools", "phangorn")

# Get a list of all installed packages
installed_packages <- installed.packages()[, "Package"]

# Check if each package is installed
are_installed <- sapply(packages_to_check, function(pkg) pkg %in% installed_packages)

packages_missing <- packages_to_check[!are_installed]

if (length(packages_missing) > 0) {
  message("Missing packages: ", paste(packages_missing, collapse = ", "))
  if ("--auto-install" %in% args) {
	  answer <- "y"
	} else {
  cat("Do you want to install all missing packages? [Y/N]: ")
  answer <- tolower(readLines("stdin", n = 1))
}
  
  if (answer == "y") {
    for (pkg in packages_missing) {
      if (pkg == "simclock") {
        if (!requireNamespace("devtools", quietly = TRUE)) {
          install.packages("devtools")
        }
        devtools::install_github("dosreislab/simclock")
      } else {
        install.packages(pkg)
      }
    }
  } else {
    message("Installation skipped.")
  }
} else {
  message("All required packages are already installed.")
}

q()


