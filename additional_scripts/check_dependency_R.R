#! /bin/env Rscript


#####################################################################
# List of packages to check
packages_to_check <- c("getopt", "seqinr", "phytools", "phangorn", "simclock")

# Get a list of all installed packages
installed_packages <- installed.packages()[, "Package"]

# Check if each package is installed
are_installed <- sapply(packages_to_check, function(pkg) pkg %in% installed_packages)

packages_missing <- packages_to_check[!are_installed]

if (length(packages_missing) > 0) {
    message("Installing missing packages: ", paste(packages_missing, collapse = ", "))
    #install.packages(packages_missing)
    for (pkg in packages_missing){
        if(pkg == 'simclock'){
            devtools::install_github("dosreislab/simclock")
        } else{
            install.packages(pkg)
        }        
    }
} else {
    q()
    #message("All packages are already installed.")
}

