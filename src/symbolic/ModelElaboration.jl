"""
Modia module for elaboration of models.

* Developer: Hilding Elmqvist, Mogram AB  
* First version: July-August 2016
* Copyright (c) 2016-2018: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)
"""
module ModelElaboration

const PrintOriginalModel = false # Needed for the moment to set partial attribute.
const PrintInstantiated = false
const PrintFlattened = true
const PrintElaborated = true
const PrintJSON = false
const fileStdOut = false

# ----------------------------------------------------------------------------

#using Debug
using ..Instantiation
using ..Execution
using Base.Meta: quot, isexpr
using Base: isapprox
using DataStructures
import ..Instantiation: GetField, This, Der, Symbolic, vars_of, Instance
import ..Execution: logTiming
using ..Synchronous
using ..BLTandPantelides
using ..BLTandPantelidesUtilities
using ..BasicStructuralTransform
using ..StructuralTransform
using ..SymbolicTransform
using ..Utilities
using ..ModiaLogging
#= Temporarily removed due to problem with PyPlot
using PyPlot
=#
import ModiaMath

using JSON

@static if VERSION < v"0.7.0-DEV.2005"
else
    using Printf
end

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

export elaborate, prettyPrint, simulateModel, simulate, skewCoords, skew, residue, residue_der, modiaCross, showExpr, transformModel, simulateMultiModeModel, allInstances
# export modiaSwitches, defineSwitch, setSwitch, showSwitches, getSwitch, 
export checkSimulation

const useIncidenceMatrix = BasicStructuralTransform.useIncidenceMatrix
const elaborate = true
#@show elaborate

const noResult = Dict{Symbol,AbstractArray{T,1} where T}()

allInstances(s) = return s

function dummyDer(e)
    f = e
    if typeof(f) == Instantiation.GetField
        f = f.name
    end

    if typeof(f) == Instantiation.Der
        f = dummyDer(f.base)
    end
    s = string("der_", f)
    return Symbol(s)
end
  
function dummyState(e)
    f = e
    if typeof(f) == Instantiation.GetField
        f = f.name
    end

    if typeof(f) == Instantiation.Der
        f = dummyState(f.base)
    end
    string(f)
end

function substituteDer(ex, dummyDerivatives, unassigned)
    if typeof(ex) == Instantiation.Der 
        newname = dummyDer(ex.base)
        dumSt = dummyState(ex)
        unass = map(string, unassigned)

        if dumSt in unass
            # println("substituteDer for unassigned: ", newname)
            ex = GetField(This(), newname)
        elseif newname in dummyDerivatives
            #  println("substituteDer: ", newname)
            ex = GetField(This(), newname)
        else 
            ex
        end

    elseif typeof(ex) == Expr 
        Expr(ex.head, [substituteDer(arg, dummyDerivatives, unassigned) for arg in ex.args]...)
    else
        ex
    end
end

function elaborateModel(flat_model)
    loglnModia("elaborate:")
    solved_model = flat_model
  
    #  print("Elaborate:          ")
    #  @time 
    (equations, variables, assign, components, Avar, Bequ, states, deriv, unassignedNames, incidenceMatrix, varSizes, varTypes, equSizes, equTypes) = StructuralTransform.transformStructurally(flat_model)
  
    if equations == nothing # Connector
        return
    end
  
    flat_model.equations = equations

    loglnModia("\nSYMBOLIC PROCESSING")
  
    if logTiming
        print("Symbolic processing:   ")
        @time solvedComponents, algebraic = differentiateSortedEquations(equations, variables, components, assign, Avar, Bequ)
    else
        solvedComponents, algebraic = differentiateSortedEquations(equations, variables, components, assign, Avar, Bequ)
    end

    if length(SymbolicTransform.dummyDerivatives.keys) > 0
        printSymbolList("\nDummy derivatives", SymbolicTransform.dummyDerivatives.keys)  
    end

    unassignedVar = []  
    solved_model.equations = []
    initSynchronousCounters()
    previousVariables = []

    for c in solvedComponents
        for (e, solved) in c
            e = substituteClocked(previousVariables, e)
            e = substituteDer(e, SymbolicTransform.dummyDerivatives.keys, unassignedVar)

            if !solved # typeof(e.args[1]) == Instantiation.Der || length(c) > 1   # 
                # Use equation
                push!(solved_model.equations, Expr(:(=), e.args[1], e.args[2]))
            else
                # Use assignment
                push!(solved_model.equations, Expr(:(:=), e.args[1], e.args[2]))
            end
        end
    end
  
    # solved_model.variables = [solved_model.variables; SymbolicTransform.dummyDerivatives]

    for d in SymbolicTransform.dummyDerivatives.keys
    # Remove der_
        v = string(d)
        if length(v) > 4 && v[1:4] == "der_"
            v = v[5:end]
        end

        if length(v) > 4 && v[1:4] == "der_"
            v = v[5:end]
        end

        vsize = varSizes[Symbol(v)]
        st = if vsize == (); 0.0 else zeros(vsize) end
        solved_model.variables[d] = Instantiation.Variable(start=st)
    end
  
    #= For test:
    for d in map(string, unassignedVar)
        @show Symbol(d)
        solved_model.variables[Symbol("der_"*d)] = Instantiation.Variable()  # add typ
        solved_model.variables[Symbol("der_der_"*d)] = Instantiation.Variable()
    end
    =#
  
    solved_model.initial_pre = []
    solved_model.F_pre = []
    solved_model.F_post = []
  
    for p in previousVariables
        push!(solved_model.initial_pre, :($(p.name) = 0))
        push!(solved_model.F_pre, :($(p.name) = 0))
        push!(solved_model.F_post, :(dummy = $(Synchronous.updatePrevious)($(p.name), $(SymbolicTransform.simulationModel_symbol), 1)))
    end
        
    if PrintElaborated 
        loglnModia("\nSTRUCTURALLY AND SYMBOLICALLY PROCESSED MODEL")
        loglnModia("  Sorted equations, := used for solved algebraic variables")
        showInstance(solved_model)
        loglnModia()
    end

    return solved_model
end


# Handle allInstances operator
function allInstances!(all, instName, inst, var, class)
    for key in keys(inst.variables) 
        if isa(inst.variables[key], Instance)
            allInstances!(all, key, inst.variables[key], var, class)
        else
            if inst.model_name == class && key == var
                push!(all, instName) 
            end
        end
    end
end

function substituteAllInstances(ex, modified_model, class, flat)
    if typeof(ex) == Expr && ex.head == :call && ex.args[1] == allInstances
        all = []
        name = string(ex.args[2].name)

@static if VERSION < v"0.7.0-DEV.2005"
        name = name[rsearchindex(name, ".") + 1:end]
else
        name = name[first(something(findlast(".", name), 0:-1)) + 1:end]
end
       
        name = Symbol(name)

        if !flat
            loglnModia("Found name: $name in $class")
            Expr(ex.head, ex.args..., class)
        else
            allInstances!(all, :(), modified_model, name, ex.args[3])
            loglnModia("all = $all")
            ref = [GetField(This(), Symbol(string(b) * "." * string(name))) for b in all]
            loglnModia("ref = $ref")
            Expr(:vect, ref...) 
        end

    elseif typeof(ex) == Expr 
        Expr(ex.head, [substituteAllInstances(arg, modified_model, class, flat) for arg in ex.args]...)
    else
        ex
    end
end

function traverseAndSubstituteAllInstances(instName, modified_model, inst, flat)
    for key in keys(inst.variables) 
        if isa(inst.variables[key], Instance)
            traverseAndSubstituteAllInstances(key, modified_model, inst.variables[key], flat)
        end
    end

    newEqu = []
    for e in inst.equations
        newE = substituteAllInstances(e, modified_model, inst.model_name, flat)
        if newE != e
            # @show e
            # @show newE
    end
        push!(newEqu, newE)
    end
    inst.equations = newEqu
end
    
function simulateModelWithOptions(model, t; options=Dict())
    opt = Dict(options)  
    ModiaLogging.setDefaultLogName(string(model.name))
    ModiaLogging.setOptions(opt)  
    StructuralTransform.setOptions(opt)  
    BasicStructuralTransform.setOptions(opt)  
    Instantiation.setOptions(opt) 
    Execution.setOptions(opt)  
  
    logSimulation = false
    if haskey(opt, :logSimulation)
        logSimulation = opt[:logSimulation]
        @show logSimulation
        delete!(opt, :logSimulation)
    end

    relTol = 1e-4
    if haskey(opt, :relTol)
        relTol = opt[:relTol]
        @show relTol
        delete!(opt, :relTol)
    end

    hev = 1e-8
    if haskey(opt, :hev)
        hev = opt[:hev]
        @show hev
        delete!(opt, :hev)
    end

    if length(keys(opt)) > 0
        println("Option(s) not found: ", keys(opt))
    end
    
    if fileStdOut
        @static if VERSION < v"0.7.0-DEV.2005"
            originalSTDOUT = STDOUT
        else
            originalSTDOUT = stdout
        end
        (outRead, outWrite) = redirect_stdout()
    end  
  
    openLogModia()
  
    start = time_ns()
    setDerAsFunction(false)
  
    if BasicStructuralTransform.logStatistics
        println("\nSimulating model: ", model.name)
    end
    loglnModia("\nSimulating model: ", model.name)
  
    if PrintOriginalModel 
        loglnModia("ORIGINAL MODEL:")
        showModel(model)
    end

    if logTiming
        print("Instantiate:           ")
        @time modified_model = instantiate(model, first(t))
    else
        modified_model = instantiate(model, first(t))
    end
    traverseAndSubstituteAllInstances(:(), modified_model, modified_model, false)

    if PrintInstantiated 
        loglnModia("INSTANTIATED MODEL:")
        showInstance(modified_model)
        loglnModia()
    end
  
    if PrintJSON
        name = string(modified_model.model_name)
        file = open(name * ".JSON", "w")
        printJSON(file, modified_model, name, name)
        close(file)
    end

    if logTiming
        print("Flatten:               ")
        @time flat_model = flatten(modified_model)
    else
        flat_model = flatten(modified_model)
    end
    traverseAndSubstituteAllInstances(:(), modified_model, flat_model, true)

    if PrintFlattened
        loglnModia("\nFLATTENED MODEL")
        showInstance(flat_model)
        loglnModia()
    end
  
    if elaborate
        solved_model = elaborateModel(flat_model)
        if solved_model == nothing
            return nothing
        end
    else
        solved_model = flat_model
    end
  
    #=
    #  JSON.json(t::This) = "this."
  
    if PrintJSONsolved
       jsonSolved = JSON.json(solved_model.equations[1])
        @show jsonSolved
    end
    =#    
    @static if VERSION < v"0.7.0-DEV.2005"
        disableSimulation = false
    else
        disableSimulation = false
    end  
  
    if !BasicStructuralTransform.newStateSelection && !disableSimulation ## Disable simulation for the moment
        loglnModia("\nSIMULATION")
        # @show incidenceMatrix
        if logTiming
            print("Code generation and simulation:         ")
            # @show solved_model t useIncidenceMatrix logSimulation
            @time res = simulate_ida(solved_model, t, if useIncidenceMatrix; incidenceMatrix else nothing end, log=logSimulation, relTol=relTol, hev=hev)
        else
            res = simulate_ida(solved_model, t, if useIncidenceMatrix; incidenceMatrix else nothing end, log=logSimulation, relTol=relTol, hev=hev)
        end
    else
        res = Dict{Symbol,AbstractArray{T,1} where T}()
    end
  
    setDerAsFunction(true)

    if fileStdOut
        #  close(outWrite)
        output = readavailable(outRead)
        close(outRead)
        redirect_stdout(originalSTDOUT)
        print(String(output))
    end
  
    if logTiming
        println(@sprintf("Total time: %0.3f", (time_ns() - start) * 1E-9), " seconds")
    end
    closeLogModia()
    @test res != noResult
    return res
end


function simulateModel(model, t; options...) 
    return simulateModelWithOptions(model, t, options=options)
end

function simulate(model, stopTime; startTime=0, options...)
    nSteps = 1000
    @static if VERSION < v"0.7.0-DEV.2005"
        t = linspace(startTime, stopTime, nSteps) 
    else
        t = range(startTime, stop=stopTime, length=nSteps)
    end
    
    return simulateModelWithOptions(model, t, options=options)
end

function print_rgb(r, g, b, t)
    println("\e[1m\e[38;2;$r;$g;$b;249m", t)
end
     
@static if VERSION < v"0.7.0-DEV.2005"
    printstyled(s; bold=false, color=:black) = print_with_color(color, s, bold=bold)
end  
"""
    function checkSimulation(mod, stopTime, observer, finalSolution; options...)
Simulates model mod until stopTime and checks that the final value of the observer variable is approximately finalSolution. The following options are available:

* `logTranslation`: determines if logging of the translation is performed
* ...
"""      
function checkSimulation(mod, stopTime, observer="", finalSolution=0.0; startTime=0, options...)
    nSteps = 1000
    @static if VERSION < v"0.7.0-DEV.2005"
        t = linspace(startTime, stopTime, nSteps) 
    else
        t = range(startTime, stop=stopTime, length=nSteps)
    end
    res = nothing
    
    try
        res = simulateModelWithOptions(mod, t, options=options)
    catch err
        st = stacktrace(catch_backtrace())
        closeLogModia()
    
        setTestStatus(false)
        println()
        println("\n----------------------\n")
        println()
            printstyled("Simulation FAILED:", bold=true, color=:red); println()
        if isa(err, ErrorException)
            printstyled(err.msg, bold=true, color=:red); println()
        elseif isa(err, UndefVarError)
            printstyled(err, bold=true, color=:red); println()
            ModiaLogging.increaseLogCategory(:(UndefinedSymbol))      
        else
            printstyled(err, bold=true, color=:red); println()
        end
        println()
        println("\n----------------------\n")
        println()
        println("Stack trace: ------------------------------------------------")
        for l in st
            println(l)
        end
        println("End stack trace: --------------- ----------------------------")
        @test false
        return noResult
    end 
    
    final = nothing
    if res != nothing && res != noResult && observer != "" && haskey(res, observer)
        obs = res[observer]
        if length(obs) > 0
            final = obs[end]
        end
    end
    
    if final != nothing
        println("final $observer = $final")
        ok = finalSolution == nothing || typeof(final) == Float64 && Base.isapprox(final, finalSolution, rtol=1.0E-3) || final == finalSolution
        @test ok
        
        if !ok
            setTestStatus(false)
            println("final solution $observer = $finalSolution")
            printstyled("Simulation NOT OK", bold=true, color=:red); println()
            println()
#= Replaced PyPlot solution with ModiaMath.plot(..) call
            figure()
            title("Simulation NOT OK in " * string(mod.name))
            plot(res["time"], res[observer])
            grid(true)
            xlabel("time [s]")
            legend([observer],  loc="center right")
=#
            ModiaMath.plot(res, observer, heading="Simulation NOT OK in " * string(mod.name))
            false
        else
            ModiaLogging.increaseLogCategory(:(CalculatedResult))   
            setTestStatus(true)
            printstyled("Simulation OK", bold=true, color=:green); println()
            println()
            true
        end
    
    else
        setTestStatus(true)
        printstyled("Simulation RAN", bold=true, color=:green); println()
        println()
        true  
        @test true
    end
    res
end

# ----------------------------


"""
Experimental code for multi-mode handling with impulses.
"""
function simulateMultiModeModel(model, t0, t1; n=1000, m=100, options...)
    opt = Dict(options)   
    ModiaLogging.setDefaultLogName(string(model.name))
    ModiaLogging.setOptions(opt)  
    StructuralTransform.setOptions(opt)  
    BasicStructuralTransform.setOptions(opt)  
    Instiantiation.setOptions(opt) 
    Execution.setOptions(opt)  
    if length(keys(opt)) >= 0
        println("Option(s) not found: ", keys(opt))
    end

    println("\nSimulating model: ", model.name)
    loglnModia("\nSimulating model: ", model.name)
  
    dt = (t1 - t0) / m
    @static if VERSION < v"0.7.0-DEV.2005"
        times1 = linspace(t0, t0 + dt, 100)  
    else
        times1 = range(t0, stop=t0 + dt, length=100)
    end
    
    model1 = flatten(instantiate(model, t0, Dict()))
    #  model1 = elaborateModel(flatten(instantiate(model, t0, Dict()))) 
    println("\nSub-simulation 1; startTime=$t0, interval=$dt")
    res1 = simulate_ida(model1, times1, if false; nothing else nothing end, log=false)
    tn = t0 + dt

    for step in 2:m
        println("\nSub-simulation $step; startTime=$tn, interval=$dt")
    
        # extract final values
 
        # finals = [name => values[end] for (name,values) in res1]
        finals = Dict(name => if length(values) > 0; values[end] else false end for (name, values) in res1)
    
        @static if VERSION < v"0.7.0-DEV.2005"
            times2 = linspace(tn, tn + dt, 100)  
        else
            times2 = range(tn, stop=tn + dt, length=100)
        end
        model2 = flatten(instantiate(model, tn, Dict())) 

        # carry over final values
        vars = vars_of(model2)
    
        for (name, final) in finals
            name = Symbol(name)
            if haskey(vars, name)
                var = vars[name]
                if isa(var, Variable)
                    var.start = final
#                   @show name final
                end
            end
        end

        res2 = simulate_ida(model2, times2, if false; nothing else nothing end, log=false)
        for (name, values) in res1
            res1[name] = vcat(res1[name], res2[name])
        end
        tn += dt
    end

    return res1
end

end