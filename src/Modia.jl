"""
Main module of Modia.

* Developer: Hilding Elmqvist, Mogram AB
* First version: December 2020
* License: MIT (expat)

"""
module Modia

const path = dirname(dirname(@__FILE__))   # Absolute path of package directory
const Version = "0.8.2"
const Date = "2022-03-08"
const modelsPath = joinpath(Modia.path, "models")

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


# Defalut switch settings
logStatistics = false
logExecution = false
logCalculations = false

useNewCodeGeneration = false
experimentalTranslation = false

using Reexport

@reexport using Unitful                 # export Unitful symbols
@reexport using DifferentialEquations   # export DifferentialEquations symbols

export ModiaBase
export CVODE_BDF, IDA
export instantiateModel, @instantiateModel, assert, stringifyDefinition
export stripUnit

export simulate!, linearize!, get_result
export @usingModiaPlot, usePlotPackage, usePreviousPlotPackage, currentPlotPackage
export resultInfo, printResultInfo, rawSignal, getPlotSignal, defaultHeading
export signalNames, timeSignalName, hasOneTimeSignal, hasSignal

export SimulationModel, measurementToString, get_lastValue
export positive, negative, previous, edge, after, reinit, pre
export initial, terminal, isInitial, isTerminal, initLinearEquationsIteration!
export get_xNames
export registerExtraSimulateKeywordArguments
export get_extraSimulateKeywordArgumentsDict

export modelToJSON, JSONToModel, writeModel, readModel

import Sundials
const  CVODE_BDF = Sundials.CVODE_BDF
const  IDA = Sundials.IDA


using Base.Meta: isexpr
using OrderedCollections: OrderedDict

using ModiaBase
using ModiaBase.Symbolic
using ModiaBase.Simplify
using ModiaBase.BLTandPantelidesUtilities
using ModiaBase.BLTandPantelides
using ModiaBase.Differentiate

import ModiaResult
import ModiaResult: usePlotPackage, usePreviousPlotPackage, currentPlotPackage
import ModiaResult: resultInfo, printResultInfo, rawSignal, getPlotSignal, defaultHeading
import ModiaResult: signalNames, timeSignalName, hasOneTimeSignal, hasSignal

import StaticArrays   # Make StaticArrays available for the tests


using  Measurements
import MonteCarloMeasurements
using JSON
#using Profile
import Test
using  TimerOutputs
using  InteractiveUtils

global to = TimerOutput()

Unitful.unit(      v::MonteCarloMeasurements.AbstractParticles{T,N}) where {T,N} = unit(T)
Unitful.upreferred(v::MonteCarloMeasurements.AbstractParticles{T,N}) where {T,N} = uconvert(upreferred(unit(v)), v)

# append! as needed in EquationAndStateInfo.jl and in CodeGeneration.jl
appendVariable!(v1::Vector{FloatType}, s::FloatType) where {FloatType} = push!(v1,s)
appendVariable!(v1::Vector{FloatType}, v2)           where {FloatType} = append!(v1,v2)


"""
    stripUnit(v)

Convert the unit of `v` to the preferred units (default are the SI units),
and then strip the unit. For details see `upreferred` and `preferunits` in
[Unitful](https://painterqubits.github.io/Unitful.jl/stable/conversion/)

The function is defined as: `stripUnit(v) = ustrip.(upreferred.(v))`.
"""
stripUnit(v) = ustrip.(upreferred.(v))


"""
     str = unitAsString( unitOfQuantity::Unitful.FreeUnits )

Return a string representation of the unit of a quantity that can be used in a unit string macro
(see also Unitful [issue 412](https://github.com/PainterQubits/Unitful.jl/issues/412)).

# Example
```
v1 = 2.0u"m/s"
v1_unit = unitAsString( unit(v1) )   # = "m*s^-1"
v2_withoutUnit = 2.0
code = :( \$v2_withoutUnit@u_str(\$v1_unit) )  # = 2.0u"m*s^-1"
v2 = eval(code)
@show v1
@show v1_unit
@show v2
```

# Notes
Transforms unit to string representation that is parseable again 
(see also Unitful [issue 412](https://github.com/PainterQubits/Unitful.jl/issues/412)).
This implementation is a hack and only works in common cases.
Implementation is performed in the following way:

1. Transform to string and display exponents on units not as Unicode superscripts (= default on macOS).
2. Replace " " by "*", since Unitful removes "*" when converting to string.
"""
unitAsString(unitOfQuantity::Unitful.FreeUnits) = replace(repr(unitOfQuantity,context = Pair(:fancy_exponent,false)), " " => "*")

"""
    quantityType = quantity(numberType, numberUnit::Unitful.FreeUnits)
    
Return the quantity type given the numberType and the numberUnit.

# Example
```julia
mutable struct Data{FloatType <: AbstractFloat}
    velocity::quantity(FloatType, u"m/s")
end

v = Data{Float64}(2.0u"mm/s")
@show v                         # v = 
```
"""
quantity(numberType, numberUnit::Unitful.FreeUnits) = Quantity{numberType, dimension(numberUnit), typeof(numberUnit)} 

quantityTypes(::Type{Unitful.Quantity{T,D,U}}) where {T,D,U} = (T,D,U)


include("EquationAndStateInfo.jl")
include("StateSelection.jl")

include("ModelCollections.jl")
include("EvaluateParameters.jl")
include("EventHandler.jl")
include("CodeGeneration.jl")
# include("GenerateGetDerivatives.jl")
include("Synchronous.jl")
include("SimulateAndPlot.jl")
include("ReverseDiffInterface.jl")
include("PathPlanning.jl")
include("JSONModel.jl")
using .JSONModel

# include("IncidencePlot.jl")
# using .IncidencePlot
# Base.Experimental.@optlevel 0

const drawIncidence = false

include("Symbolic.jl")
include("ModiaLang.jl")

end