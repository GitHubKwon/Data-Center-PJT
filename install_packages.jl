# Run once, before anything else:  julia install_packages.jl
# Installs the Julia packages the project uses (takes a few minutes the first time).
using Pkg
Pkg.add(["JuMP", "HiGHS", "CSV", "DataFrames", "XLSX", "Plots"])
println("Packages installed. Next: julia toy_example.jl")
