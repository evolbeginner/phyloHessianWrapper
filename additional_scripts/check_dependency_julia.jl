#! /bin/env julia

############################################################
using Pkg

############################################################
# Function to check if a package is installed
function check_package_installed(pkg_name)
    installed_packages = Vector()
    redirect_stderr(devnull) do
        installed_packages = keys(Pkg.installed())
    end
    return pkg_name in installed_packages
end

############################################################
try
    @assert VERSION >= v"1.0.0-dev" "Julia version must be greater than 1.0.0-dev"
catch e
    if e isa AssertionError
        println("Assertion failed: ", e)
    else
        rethrow(e)
    end
end

############################################################
# List of package names to check
package_names = [
    "ArgParse",
    "Distributions",
    #"ConcurrentCollections",
    #"CSV",
    #"DataFrames",
    "DelimitedFiles",
    #"FiniteDiff",
    #"Folds",
    #"StaticArrays",
    "QuadGK",
    "StatsFuns"
]

# Check each package
pkgs_not_found = String[]
for pkg in package_names
    if !check_package_installed(pkg)
        push!(pkgs_not_found, pkg)
    end
end

auto_install = "--auto-install" in ARGS

if !isempty(pkgs_not_found)
    println("Julia package(s) $pkgs_not_found not found!")
    if auto_install
        println("Auto-installing missing packages...")
        Pkg.add(pkgs_not_found)
    else
        print("Do you wanna install all? [Y/N]: ")
        answer = readline()
        if lowercase(answer) == "y"
            println("Installing missing packages...")
            Pkg.add(pkgs_not_found)
        else
            println("Installation skipped.")
        end
    end
end

