"""
Main module of Modia.

* Developer: Hilding Elmqvist, Mogram AB and Martin Otter, DLR
* First version: December 2020
* License: MIT (expat)

"""
module Modia

const path = dirname(dirname(@__FILE__))   # Absolute path of package directory
const Version = "0.12.0"
const Date = "2023-06-04"
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

export InstantiatedModel, measurementToString, get_lastValue, getLastValue, getStateNames
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
    @strippedPositive!(path, name)

Convert `name` to its preferred units (default are the SI units), strip units and check that value is positive.
In case of error, use `string(name)` and `path` in the error message:

# Example
```
using Unitful
L1 =  2.0u"m"
@strippedPositive!("insulatedRod", L1)
# L1 = 2.0

L2 = -2.0u"m"
@strippedPositive!("insulatedRod", L2)
# error message:
# Error from
#   insulatedRod = ...(..., L2 = -2.0u"m",...): L2 > 0 required.)
```
"""
macro strippedPositive!(path, name)
    nameAsString = string(name)
    expr = :( $name = strippedPositive!($path, $(esc(name)), $nameAsString) )   
    return expr 
end      
strippedPositive!(path::String, value, name) = stripUnit(value) > 0 ? stripUnit(value) : error("\nError from\n   $path = ...(..., $name = $value, ...): $name > 0 required")


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
include("Statistics.jl")
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

import SnoopPrecompile

SnoopPrecompile.@precompile_all_calls begin
    #= If two models are used, this gives a warning "incremental compilation broken" -> precompile only one model.
        FirstOrder = Model(
            T = 0.2u"s",
            x = Var(init=0.3),
            equations = :[u = sin(time/u"s"),
                        T * der(x) + x = u,
                        y = 2*x]
        )
        firstOrder = @instantiateModel(FirstOrder,logFile=false)
        simulate!(firstOrder, Tsit5(), stopTime = 1.0, merge = Map(T = 0.4u"s", x = 0.9))
    =#

        TwoInertiasAndIdealGear = Model(
            J1 = 0.0025,
            J2 = 170,
            r  = 105,
            tau_max = 1,
            phi1 = Var(init = nothing),
            w1   = Var(init = nothing),
            phi2 = Var(init = 0.5),
            w2   = Var(init = 0.0),
            equations = :[
                tau = sin(time/u"s"),
                w1 = der(phi1),
                J1*der(w1) = tau - tau1,
                phi1   = r*phi2,
                r*tau1 = tau2,
                w2 = der(phi2),
                J2*der(w2) = tau2
            ]
        )
        twoInertiasAndIdealGear = @instantiateModel(TwoInertiasAndIdealGear, unitless=true, logFile=false)
        simulate!(twoInertiasAndIdealGear, stopTime = 4.0, useRecursiveFactorizationUptoSize = 500)
end

end