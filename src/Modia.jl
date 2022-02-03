#=
Modia main module.

* Developers: Hilding Elmqvist (Mogram) and Martin Otter (DLR-SR)
* Copyright (c) 2016-2021: Hilding Elmqvist, Martin Otter
* License: MIT (expat)
=#

module Modia

const path = dirname(dirname(@__FILE__))   # Absolute path of package directory
const Version = "0.6.0"
const Date = "2022-02-03"

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
@reexport using ModiaLang.Unitful
@reexport using ModiaLang.DifferentialEquations
@reexport using Modia3D


const modelsPath = joinpath(ModiaLang.path, "models")


end
