"""
Main module of Modia.

* Developer: Hilding Elmqvist, Mogram AB
* First version: December 2020
* License: MIT (expat)

"""
module Modia

const path = dirname(dirname(@__FILE__))   # Absolute path of package directory
const Version = "0.9.2"
const Date = "2023-05-30"
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

@reexport using SignalTables            # export SignalTables symbols
@reexport using Unitful                 # export Unitful symbols
@reexport using DifferentialEquations   # export DifferentialEquations symbols
import SignalTables: AvailablePlotPackages

"""
    @usingModiaPlot()
    
Execute `using XXX`, where `XXX` is the Plot package that was activated with `usePlotPackage(plotPackage)`.
So this is similar to @usingPlotPackage (from SignalTables, that is reexported from Modia).

There is, however, a difference when XXX = "SilentNoPlot": 

- @usingPlotPackage() executes `using SignalTables.SilentNoPlot` and therefore requires that package `SignalTables` is available in your environment.
- @usingModiaPlot() executes `using Modia.SignalTables.SilentNoPlot` and therefore requires that package `Modia` is available in your environment.

Therefore, when working with Modia it is better to use `@usingModiaPlot()`.
"""
macro usingModiaPlot()
    if haskey(ENV, "SignalTablesPlotPackage")
        PlotPackage = ENV["SignalTablesPlotPackage"]
        if !(PlotPackage in AvailablePlotPackages)
            @info "ENV[\"SignalTablesPlotPackage\"] = \"$PlotPackage\" is not supported!. Using \"SilentNoPlot\"."
            @goto USE_NO_PLOT
        elseif PlotPackage == "NoPlot"
            @goto USE_NO_PLOT
        elseif PlotPackage == "SilentNoPlot"
            expr = :( using Modia.SignalTables.SilentNoPlot )
            return esc( expr )
        else
            PlotPackage = Symbol("SignalTablesInterface_" * PlotPackage)
            expr = :(using $PlotPackage)
            println("$expr")
            return esc( :(using $PlotPackage) )
        end

    elseif haskey(ENV, "MODIA_PLOT_PACKAGE")
        PlotPackage = ENV["MODIA_PLOT_PACKAGE"]
        if !(PlotPackage in AvailablePlotPackages)
            @info "ENV[\"MODIA_PLOT_PACKAGE\"] = \"$PlotPackage\" is not supported!. Using \"SilentNoPlot\"."
            @goto USE_NO_PLOT
        elseif PlotPackage == "NoPlot"
            @goto USE_NO_PLOT
        elseif PlotPackage == "SilentNoPlot"
            expr = :( using Modia.SignalTables.SilentNoPlot )
            return esc( expr )
        else
            PlotPackage = Symbol("SignalTablesInterface_" * PlotPackage)
            expr = :(using $PlotPackage)
            println("$expr")
            return esc( :(using $PlotPackage) )
        end

    else
        @info "No plot package activated. Using \"SilentNoPlot\"."
        @goto USE_NO_PLOT
    end

    @label USE_NO_PLOT
    expr = :( using Modia.SignalTables.SilentNoPlot )
    println("$expr")
    return esc( expr )
end
export @usingModiaPlot

export ModiaBase
export CVODE_BDF, IDA
export instantiateModel, @instantiateModel, assert, stringifyDefinition
export stripUnit

export simulate!, linearize!, get_result
export hasParameter, getParameter, getEvaluatedParameter
export showParameters, showEvaluatedParameters

export SimulationModel, measurementToString, get_lastValue, getLastValue, getStateNames
export positive, negative, previous, edge, after, reinit, pre
export initial, terminal, isInitial, isTerminal, initLinearEquationsIteration!
export get_xNames
export registerExtraSimulateKeywordArguments
export get_extraSimulateKeywordArgumentsDict

export modelToJSON, JSONToModel, writeModel, readModel

import Sundials
const  CVODE_BDF = Sundials.CVODE_BDF
const  IDA = Sundials.IDA


# Deprecated functions - only provided for backwards compatibility
export signalNames, timeSignalName, hasOneTimeSignal, printResultInfo


using Base.Meta: isexpr
using OrderedCollections: OrderedDict

using ModiaBase
using ModiaBase.Symbolic
using ModiaBase.Simplify
using ModiaBase.BLTandPantelidesUtilities
using ModiaBase.BLTandPantelides
using ModiaBase.Differentiate

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
    str = modelPathAsString(modelPath::Union{Expr,Symbol,Nothing})

Return modelPath of submodel as string.
"""
modelPathAsString(modelPath::Union{Expr,Symbol,Nothing}) = isnothing(modelPath) ? "" : string(modelPath)

include("EquationAndStateInfo.jl")
include("Result.jl")
include("StateSelection.jl")

include("ModelCollections.jl")
include("EventHandler.jl")
include("CodeGeneration.jl")
include("EvaluateParameters.jl")

# include("GenerateGetDerivatives.jl")
include("Synchronous.jl")
include("SimulateAndPlot.jl")
include("SignalTablesInterface.jl")
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