#=
Modia main module.

* Developers: Hilding Elmqvist (Mogram) and Martin Otter (DLR-SR)
* Copyright (c) 2016-2021: Hilding Elmqvist, Martin Otter
* License: MIT (expat)
=#

module Modia

const path = dirname(dirname(@__FILE__))   # Absolute path of package directory

const Version = "0.7.0"
const Date = "2021-04-25"

#println(" \n\nWelcome to Modia - Dynamic MODeling and Simulation in julIA")
print(" \n\nWelcome to ")
print("Mod")
printstyled("ia", bold=true, color=:red)
print(" - ")
printstyled("Dynamic ", color=:light_black)
print("Mod")
printstyled("eling and Simulation with Jul", color=:light_black)
printstyled("ia", bold=true, color=:red)

println()
println("Version $Version ($Date)")

using Reexport

@reexport using ModiaLang
@reexport using Unitful
@reexport using DifferentialEquations

include("../models/AllModels.jl")

end
