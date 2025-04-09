#! /bin/env julia


############################################################
using Pkg

# Function to check if a package is installed
function check_package_installed(pkg_name)
    installed_packages = Pkg.TOML.parsefile(Pkg.project().path)["deps"]
    return haskey(installed_packages, pkg_name)
end

# List of package names to check
package_names = [
    "ArgParse",
    "ConcurrentCollections",
    "CSV",
    "DataFrames",
    "DelimitedFiles",
    "FiniteDiff",
    "Folds",
    "StaticArrays",
    "StatsFuns"
]

# Check each package
for pkg in package_names
    if ! check_package_installed(pkg)
        println("Julia package $pkg not found! Try 'using Pkg; Pkg.add($pkg)' to install it.")
    end
end



