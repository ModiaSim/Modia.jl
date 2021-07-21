"""
Main module of TinyModia.

* Developer: Hilding Elmqvist, Mogram AB
* First version: December 2020
* License: MIT (expat)

"""
module TinyModia

# Defalut switch settings
logStatistics = false
logExecution = false
logCalculations = false

useNewCodeGeneration = false

using Reexport

export instantiateModel, @instantiateModel, assert, stringifyDefinition

using Base.Meta: isexpr
using OrderedCollections: OrderedDict

using ModiaBase.Symbolic
using ModiaBase.Simplify
using ModiaBase.BLTandPantelidesUtilities
using ModiaBase.BLTandPantelides
using ModiaBase.Differentiate
using ModiaBase
export ModiaBase

@reexport using Unitful
using  Measurements
import MonteCarloMeasurements
using JSON
#using Profile
using TimerOutputs
using InteractiveUtils

global to = TimerOutput()

export stripUnit
"""
    stripUnit(v)
    
Convert the unit of `v` to the preferred units (default are the SI units),
and then strip the unit. For details see `upreferred` and `preferunits` in 
[Unitful](https://painterqubits.github.io/Unitful.jl/stable/conversion/)

The function is defined as: `stripUnit(v) = ustrip(upreferred.(v))`.
"""
stripUnit(v) = ustrip.(upreferred.(v))


include("ModelCollections.jl")
include("EvaluateParameters.jl")
include("EventHandler.jl")
include("CodeGeneration.jl")
# include("GenerateGetDerivatives.jl")
include("Synchronous.jl")
include("SimulateAndPlot.jl")
include("ReverseDiffInterface.jl")
include("PathPlanning.jl")

# include("IncidencePlot.jl")
# using .IncidencePlot
# Base.Experimental.@optlevel 0

const drawIncidence = false

const path = dirname(dirname(@__FILE__))   # Absolute path of package directory

const Version = "0.8.0-dev"
const Date = "2021-07-21"

#println(" \n\nWelcome to Modia - Dynamic MODeling and Simulation in julIA")
print(" \n\nWelcome to ")
printstyled("Tiny", color=:light_black)
print("Mod")
printstyled("ia", bold=true, color=:red)
print(" - ")
printstyled("Dynamic ", color=:light_black)
print("Mod")
printstyled("eling and Simulation with Jul", color=:light_black)
printstyled("ia", bold=true, color=:red)

println()
println("Version $Version ($Date)")

# ----------------------------------------------------------------------------------------------

function assert(condition, message)
    if ! condition
        println("ASSERTION: ", message)
    end
end

function printArray(array, heading; order=1:length(array), log=false, number=true, indent=false, spacing=" ")
    if length(array) > 0 && log
        println(heading)
        n = 0
        for v in array
            n += 1
            if number
                print(lpad(order[n], 4), ": ")
            end
            if indent
                print(("|"*spacing)^(n-1))
            end
            println(replace(replace(string(v), "\n" => " "), "  " => " "))
        end
    end
end

function performConsistencyCheck(G, Avar, vActive, parameters, unknowns, states, equations, log=false)
    nUnknowns = length(unknowns) - length(states)
    nEquations = length(equations)

    if nEquations > 0 && nUnknowns != nEquations
        printstyled("The number of unknowns ($nUnknowns) and equations ($nEquations) is not equal!\n", bold=true, color=:red)
        printArray(parameters, "Parameters:", log=true)
        printArray(states, "Potential states:", log=true)
        printArray(setdiff(unknowns,states), "Unknowns:", log=true)
        printArray(equations, "Equations:", log=true)
    end

    ok = nUnknowns == nEquations

    # Check consistency
    EG::Array{Array{Int,1},1} = [G; buildExtendedSystem(Avar)]
    assignEG = matching(EG, length(Avar))
#        @assert all(assignEG .> 0) "  The DAE is structurally singular\n  $assignEG."
    assigned = [a > 0 for a in assignEG]
    if ! (all(assigned)) || ! ok
        if nUnknowns > nEquations
            missingEquations = nUnknowns - nEquations
            printstyled("\nThe DAE has no unique solution, because $missingEquations equation(s) missing!\n", bold=true, color=:red)
        elseif nEquations > nUnknowns
            tooManyEquations = nEquations - nUnknowns
            printstyled("\nThe DAE has no unique solution, because $tooManyEquations equation(s) too many\n", bold=true, color=:red)
        else
            printstyled("\nThe DAE has no unique solution, because it is structurally inconsistent\n" *
                        "  (the number of equations is identical to the number of unknowns)\n", bold=true, color=:red)
        end 
        
        singular = true
        hequ = []
        # Create dummy equations for printing: integrate(x, der(x)) = 0
        for i in 1:length(Avar)
            a = Avar[i]
            if a > 0
                n = unknowns[i]
                der = unknowns[a]
                push!(hequ, "integrate($n, $der) = 0")
            end
        end
        equationsEG = vcat(equations, hequ)
        println("\nAssigned equations after consistency check:")
        printAssignedEquations(equationsEG, unknowns, 1:length(EG), assignEG, Avar, fill(0, length(EG)))
        unassigned = []
        for i in eachindex(assignEG)
            if assignEG[i] == 0
                push!(unassigned, unknowns[i])
            end
        end
        printArray(unassigned, "Not assigned variables after consistency check:", log=true)
        if nUnknowns > nEquations
            missingEquations = nUnknowns - nEquations
            error("Aborting, because $missingEquations equation(s) missing!")
        elseif nEquations > nUnknowns
            tooManyEquations = nEquations - nUnknowns
            error("Aborting, because $tooManyEquations equation(s) too many!")   
        else
            error("Aborting, because DAE is structurally inconsistent")
        end                
    elseif log
        println("  The DAE is structurally nonsingular.")
    end
end



function assignAndBLT(equations, unknowns, parameters, Avar, G, states, logDetails, log=false, logTiming=false)
    vActive = [a == 0 for a in Avar]

    if drawIncidence
        showIncidence(equations, unknowns, G, vActive, [], [], [], [[e] for e in 1:length(G)])
    end

    if false # logAssignmentBeforePantelides
        @show G length(Avar) vActive
        assign = matching(G, length(Avar), vActive)
        @show assign

        println("\nAssigned equations:")
        Bequ = fill(0, length(G))
        printAssignedEquations(equations, unknowns, 1:length(G), assign, Avar, Bequ)

        printUnassigned(equations, unknowns, assign, Avar, Bequ, vActive)
    end

    if logTiming
        println("\nCheck consistency of equations by matching extended equation set")
        @timeit to "performConsistencyCheck" performConsistencyCheck(G, Avar, vActive, parameters, unknowns, states, equations, log)
    else
        performConsistencyCheck(G, Avar, vActive, parameters, unknowns, states, equations, log)
    end

    if logTiming
        println("\nIndex reduction (Pantelides)")
        @timeit to "pantelides!" assign, Avar, Bequ = pantelides!(G, length(Avar), Avar)
    else
        assign, Avar, Bequ = pantelides!(G, length(Avar), Avar)
    end
    if logDetails
        @show assign Avar Bequ
    end

    moreVariables = length(Avar) > length(unknowns)
    unknowns = vcat(unknowns, fill(:x, length(Avar) - length(unknowns)))
    for i in 1:length(Avar)
        if Avar[i] > 0
            unknowns[Avar[i]] = derivative(unknowns[i], keys(parameters))
        end
    end

    moreEquations = length(Bequ) > length(equations)
    equations = vcat(equations, fill(:(a=b), length(Bequ) - length(equations)))
    for i in 1:length(Bequ)
        if Bequ[i] > 0
            if log
                println("Differentiating: ", equations[i])
            end
            equations[Bequ[i]] = simplify(derivative(equations[i], keys(parameters)))
            if log
                println("  Differentiated equation: ", equations[Bequ[i]])
            end
        end
    end

    if moreVariables && log
        printArray(unknowns, "Unknowns after index reduction:", log=logDetails)
    end
    if moreEquations && log
        printArray(equations, "Equations after index reduction:", log=logDetails)
    end

    HG = Array{Array{Int64,1},1}()
    originalEquationIndex = []
    highestDerivativeEquations = []
    for i in 1:length(G)
        if Bequ[i] == 0
            push!(HG, G[i])
            push!(highestDerivativeEquations, equations[i])
            push!(originalEquationIndex, i)
        end
    end

    vActive = [a == 0 for a in Avar]
    if logTiming
        println("Assign")
        @timeit to "matching" assign = matching(HG, length(Avar), vActive)
    else
        assign = matching(HG, length(Avar), vActive)
    end
    if logTiming
        println("BLT")
        @timeit to "BLT" bltComponents = BLT(HG, assign)
    else
       bltComponents = BLT(HG, assign)
    end
    if logDetails
        @show HG
        @show bltComponents
    end

    if log
        println("\nSorted highest derivative equations:")
        printSortedEquations(highestDerivativeEquations, unknowns, bltComponents, assign, Avar, Bequ)
    end

    if drawIncidence
        assignedVar, unAssigned = invertAssign(assign)
        showIncidence(equations, unknowns, HG, vActive, [], assign, assignedVar, bltComponents, sortVariables=true)
    end

    for b in bltComponents
        if length(b) > 1 && logDetails
            println("\nTearing of BLT block")
            @show b
            es::Array{Int64,1} = b
            assignedVar, unAssigned = invertAssign(assign)
            vs = assignedVar[b]
            td = TearingSetup(HG, length(unknowns))
            (eSolved, vSolved, eResidue, vTear) = tearEquations!(td, (e,v) -> v in HG[e], es, vs)
            @show vTear eResidue vSolved eSolved
            printArray(unknowns[vTear], "torn variables:", log=logDetails)
            printArray(equations[eResidue], "residue equations:", log=logDetails)
            printArray(unknowns[vSolved], "solved variables:", log=logDetails)
            printArray(equations[eSolved], "solved equations:", log=logDetails)
        end
    end

    blt = Vector{Int}[]  # TODO: Change type
    for c in bltComponents
        push!(blt, originalEquationIndex[c])
    end
    for i in eachindex(assign)
        if assign[i] != 0
            assign[i] = originalEquationIndex[assign[i]]
        end
    end
    return (unknowns, equations, G, Avar, Bequ, assign, blt, parameters)
end

function setAvar(unknowns)
    variablesIndices = OrderedDict(key => k for (k, key) in enumerate(unknowns))
    Avar = fill(0, length(unknowns))
    states = []
    derivatives = []
    for i in 1:length(unknowns)
        v = unknowns[i]
        if isexpr(v, :call)
            push!(derivatives, v)
            x = v.args[2]
            ix = variablesIndices[x] 
            Avar[ix] = i
            push!(states, unknowns[ix])
        end
    end
    return Avar, states, derivatives
end

mutable struct ModelStructure
    parameters::OrderedDict{Any,Any}
    mappedParameters::OrderedDict{Any,Any} 
    init::OrderedDict{Any,Any}
    start::OrderedDict{Any,Any}
    variables::OrderedDict{Any,Any}
    flows::OrderedDict{Any,Any}
    inputs::OrderedDict{Any,Any}
    outputs::OrderedDict{Any,Any}
#    components
#    extends
    equations::Array{Expr,1}
end

ModelStructure() = ModelStructure(OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), Expr[])

function printDict(label, d)
    if length(d) > 0
        print(label)
        for (k,v) in d
            print(k)
            if v != nothing
                print(" => ", v)
            end
            print(", ")
        end
        println()
    end
end

function printModelStructure(m; label="", log=false)
    if log
        println("ModelStructure ", label)
        printDict("  parameters:  ", m.parameters)
        printDict("  mappedParameters:  ", m.mappedParameters)
        printDict("  init:        ", m.init)
        printDict("  start:       ", m.start)
        printDict("  variables:   ", m.variables)
        printDict("  flows:       ", m.flows)
        printDict("  inputs:      ", m.inputs)
        printDict("  outputs:     ", m.outputs)
        println("  equations: ")
        for e in m.equations
            println("    ", e)
        end
    end
end

prependDict(dict, prefix) = OrderedDict([prepend(k, prefix) => prepend(v, prefix) for (k,v) in dict])


function mergeModelStructures(parent::ModelStructure, child::ModelStructure, prefix)
    merge!(parent.parameters, child.parameters)
    parent.mappedParameters[prefix] = child.mappedParameters

    merge!(parent.init, child.init)
    parent.mappedParameters[prefix] = child.mappedParameters

    merge!(parent.start, child.start)
    parent.mappedParameters[prefix] = child.mappedParameters

    merge!(parent.variables, child.variables)
    merge!(parent.flows, child.flows)

    push!(parent.equations, prepend(child.equations, prefix)...)
end

function replaceLinearIntegerEquations(Gint, eInt, GcInt, unknowns, equations)
    for i = 1:length(Gint)
        e = Gint[i]
        if length(e) > 0
            # Construct equation
            rhs = 0
            for j = 1:length(e)
                v  = e[j]
                vc = GcInt[i][j]
                rhs = add(rhs, mult(vc, unknowns[v]))
            end
            equ = :(0 = $rhs)
            equations[eInt[i]] = equ
        else
            equations[eInt[i]] = :(0 = 0)
        end
    end
end

function performAliasReduction(unknowns, equations, Avar, logDetails, log)

    @timeit to "enumerate(unknowns)" variablesIndices = OrderedDict(key => k for (k, key) in enumerate(unknowns))

    linearVariables = []
    nonlinearVariables = []
    linearEquations = Vector{Int}()
    nonlinearEquations = []
    G = Vector{Int}[] # Array{Array{Int64,1},1}()
    Gint = Array{Array{Int64,1},1}()
    GcInt = Array{Array{Int64,1},1}()

    @timeit to "build graph information" for i in 1:length(equations)
        e = equations[i]
        nameIncidence, coefficients, rest, linear = getCoefficients(e)
        incidence = [variablesIndices[n] for n in nameIncidence if n in keys(variablesIndices)]
        unique!(incidence)
        push!(G, incidence)

        if linear && all([c == 1 || c == -1 for c in coefficients]) && rest == 0 && all(n in keys(variablesIndices) for n in nameIncidence)
            push!(linearVariables, nameIncidence...)
            push!(linearEquations, i)
            push!(Gint, incidence)
            push!(GcInt, coefficients)
        else
            push!(nonlinearVariables, nameIncidence...)
            push!(nonlinearEquations, e)
        end
    end

    if logDetails
        @show G
        @show Avar
    end

    @timeit to "unique!(linearVariables)" unique!(linearVariables)
    @timeit to "unique!(nonlinearVariables)" unique!(nonlinearVariables)
    @timeit to "setdiff" vSolveInLinearEquations = setdiff(linearVariables, nonlinearVariables)
    if logDetails
        printArray(linearVariables, "linearVariables:", number=false)
        printArray(equations[linearEquations], "linearEquations:", number=false)
#        printArray(nonlinearVariables, "nonlinearVariables:", number=false)
#        printArray(nonlinearEquations, "nonlinearEquations:", number=false)
        @show linearEquations vSolveInLinearEquations Gint GcInt
    end

    GCopy = deepcopy(G)
    AvarCopy = deepcopy(Avar)
    @timeit to "simplifyLinearIntegerEquations!" (vEliminated, vProperty, nvArbitrary, redundantEquations) = simplifyLinearIntegerEquations!(GCopy, linearEquations, GcInt, AvarCopy)

    if logDetails
        @show vEliminated vProperty nvArbitrary redundantEquations
    end
    if log
        printSimplifiedLinearIntegerEquations(GCopy, linearEquations, GcInt, vEliminated, vProperty, nvArbitrary, redundantEquations, v->string(unknowns[v]), printTest=false)
    end

    substitutions = OrderedDict()
    @timeit to "build substitutions" if length(vEliminated) - nvArbitrary > 0
        for v in vEliminated[nvArbitrary+1:end]
            if isZero(vProperty, v)
                substitutions[unknowns[v]] = 0.0
            elseif isAlias(vProperty, v)
                substitutions[unknowns[v]] = unknowns[alias(vProperty,v)]
            else
                substitutions[unknowns[v]] = :(- $(unknowns[negAlias(vProperty,v)]))
            end
        end
    end

    @timeit to "replaceLinearIntegerEquations" replaceLinearIntegerEquations(GCopy[linearEquations], linearEquations, GcInt, unknowns, equations)
#    printArray(equations, "new equations:", log=log)

    reducedEquations = []
    @timeit to "substitute" for ei in 1:length(equations)
        e = equations[ei]
        if e != :(0 = 0)
            push!(reducedEquations, substitute(substitutions, e))
        end
    end

    reducedUnknowns = []
    for vi in 1:length(unknowns)
        if isNotEliminated(vProperty, vi)
            push!(reducedUnknowns, unknowns[vi])
        end
    end

    reducedVariablesIndices = OrderedDict(key => k for (k, key) in enumerate(reducedUnknowns))
    reducedG = Vector{Int}[] # Array{Array{Int64,1},1}()
    @timeit to "build reducedG" for i in 1:length(reducedEquations)
        e = reducedEquations[i]
        nameIncidence, coefficients, rest, linear = getCoefficients(e)
        incidence = [reducedVariablesIndices[n] for n in nameIncidence if n in keys(reducedVariablesIndices)]
        unique!(incidence)
        push!(reducedG, incidence)
    end

    reducedAvar, reducedStates, reducedDerivatives = setAvar(reducedUnknowns)
    return reducedEquations, reducedUnknowns, reducedAvar, reducedG, reducedStates, vEliminated, vProperty
end


function stateSelectionAndCodeGeneration(modStructure, name, modelModule, FloatType, init, start, inputs, outputs, vEliminated, vProperty, unknownsWithEliminated, mappedParameters;
    unitless=false, logStateSelection=false, logCode=false, logExecution=false, logCalculations=false, logTiming=false, evaluateParameters=false)
    (unknowns, equations, G, Avar, Bequ, assign, blt, parameters) = modStructure

    function getSolvedEquationAST(e, v)
        (solution, solved) = solveEquation(equations[e], unknowns[v])
        unknown = solution.args[1]
        expr = solution.args[2]
###        solution = :(_variables[$(v+1)] = $unknown = $expr)
        equ = string(equations[e])
        @assert solved "Equation not solved for $(unknowns[v]): $(equations[e])"
#        sol = string(solution)
#        solution = prepend(makeDerivativeVar(solution, components), :m)
        solution = substituteForEvents(solution)
        if false
    #        solution = :(begin println("Executing: ", $sol); $solution end)
            if sol != equ
                solution = :(try $solution; catch e; println("\n\nFailure executing: ", $sol); println("Original equation: ", $equ, "\n\n"); rethrow(e) end)
            else
                solution = :(try $solution; catch e; println("\n\nFailure executing: ", $sol, "\n\n"); rethrow(e) end)
            end
    #        solution = :(try $solution; catch e; println("Failure executing: ", $sol); println(sprint(showerror,e)); rethrow(e) end)
    #        solution = :(try $solution; catch e; println("Failure executing: ", $sol); printstyled(stderr,"ERROR: ", bold=true, color=:red);
    #        printstyled(stderr,sprint(showerror,e), color=:light_red); println(stderr); end)
        end
        solution = makeDerVar(solution, parameters, inputs, evaluateParameters)
        if logCalculations
            var = string(unknowns[v])
            solutionString = string(solution)
            return :(println("Calculating: ", $solutionString); $solution; println("  Result: ", $var, " = ", upreferred.($(solution.args[1]))))
        else
            return solution
        end
    end

    function getResidualEquationAST(e, residualName)
        eq = equations[e] # prepend(makeDerivativeVar(equations[e], components), :m)
        eqs = sub(eq.args[2], eq.args[1])
        resid = makeDerVar(:(ustrip($eqs)), parameters, inputs, evaluateParameters)
        residual = :($residualName = $resid)
        residString = string(resid)
        if logCalculations
            return :(println("Calculating residual: ", $residString); $residualName = $resid; println("  Residual: ", $residualName) )
#            return makeDerVar(:(dump($(makeDerVar(eq.args[2]))); dump($(makeDerVar(eq.args[1]))); $residual; println($residualName, " = ", upreferred.(($(eq.args[2]) - $(eq.args[1]))))))
        else
            return residual
        end
    end

    function isLinearEquation(e_original, v_original)
        equ = equations[e_original]
        var = unknowns[v_original]
        return isLinear(equ, var)
    end

    function isSolvableEquation(e_original, v_original)
        equ = equations[e_original]
        var = unknowns[v_original]
        (solution, solved) = solveEquation(equ, var)
        return solved
    end

    hasParticles(value) = typeof(value) <: MonteCarloMeasurements.StaticParticles ||
                          typeof(value) <: MonteCarloMeasurements.Particles

    function var_unit(v)
        var = unknowns[v]
        if var in keys(init)
            value = eval(init[var])
        elseif var in keys(start)
            value = eval(start[var])
        else
            value = 0.0
        end
        if hasParticles(value)  # Units not yet support for particles
            return ""
        end
        # if length(value) == 1
        if ! (typeof(value) <: Array)        
            un = unit(value)
        else
            un = unit.(value)   # un = [unit(v) for v in value]  # unit.(value) does not work for MonteCarloMeasurements
            @assert all([un[i] == un[1] for i in 2:length(un)]) "The unit of all elements of state vector must be equal: $var::$(value)"
            un = un[1]
        end
        return replace(string(un), " " => "*") # Fix since Unitful removes * in unit strings
    end

    function var_has_start(v_original)
        return unknowns[v_original] in keys(init) || unknowns[v_original] in keys(start)
    end

    function var_fixed(v_original)
        return unknowns[v_original] in keys(init)
    end

    function var_length(v_original)
        v = unknowns[v_original]
        if v in keys(init)
            value = eval(init[v])
            len = hasParticles(value) ? 1 : length(value)
        elseif v in keys(start)
            value = eval(start[v])
            len = hasParticles(value) ? 1 : length(value)            
        else
            len = 1
        end
        return len
    end

    function eq_equation(e)
        return replace(replace(string(equations[e]), "\n" => " "), "  " => " ")
    end

    # ----------------------------------------------------------------------------

    solvedAST = []
    Gsolvable = copy(G)
    juliaVariables  = [Symbol(u) for u in unknowns]
    stringVariables = [String(Symbol(v)) for v in unknowns]
    allEquations    = equations
    #    @show stringVariables equations G blt assign Avar Bequ Gsolvable


    stateSelectionFunctions = StateSelectionFunctions(
        var_name               = v -> stringVariables[v],
        var_julia_name         = v -> juliaVariables[v],
        var_unit               = var_unit,
        var_has_start          = var_has_start,
        var_fixed              = var_fixed,
        var_nominal            = v_original -> NaN,
        var_unbounded          = v_original -> false,
        var_length             = var_length,
        equation               = eq_equation,
        isSolvableEquation     = isSolvableEquation,
        isLinearEquation       = isLinearEquation,
        getSolvedEquationAST   = getSolvedEquationAST,
        getResidualEquationAST = getResidualEquationAST,
        showMessage            = (message;severity=0,from="???",details="",variables=Int[],equations=Int[]) ->
                                    ModiaBase.showMessage(message, severity, from, details,
                                                          stringVariables[variables],
                                                          string.(allEquations[equations]))
        )

    if length(equations) == 0
        return nothing
    end

    if logTiming
        println("Get sorted and solved AST")
        @timeit to "getSortedAndSolvedAST" (success,AST,equationInfo) = getSortedAndSolvedAST(
            G, Array{Array{Int64,1},1}(blt), assign, Avar, Bequ, stateSelectionFunctions;
            unitless, log = logStateSelection, modelName = name)
    else
        (success,AST,equationInfo) = getSortedAndSolvedAST(
            G, Array{Array{Int64,1},1}(blt), assign, Avar, Bequ, stateSelectionFunctions;
            unitless, log = logStateSelection, modelName = name)
    end

    if ! success
        error("Aborting")
    end

    if useNewCodeGeneration
        AST, newFunctions = prepareGenerateGetDerivatives(AST, modelModule, functionSize, juliaVariables)
    else
        newFunctions = Vector{Expr}()
    end

#    @show AST
#    @show equationInfo

    selectedStates = [v.x_name_julia for v in equationInfo.x_info]

    startValues = []
    for v in selectedStates
        v = Meta.parse(string(v))
        if v in keys(init)
            value = eval(init[v])
#            value = [0.0 0.0]  # TODO: Fix
#            @show value
            if hasParticles(value)
                push!(startValues, value)
            else
                push!(startValues, value...)
            end
        elseif v in keys(start)
            value = eval(start[v])
            if hasParticles(value)
                push!(startValues, value)
            else
                push!(startValues, value...)
            end
        else
            push!(startValues, 0.0)
        end
    end
    if logCode
        println("startValues = ", startValues)
    end

    vSolvedWithInit = equationInfo.vSolvedWithFixedTrue
    vSolvedWithInitValuesAndUnit = OrderedDict{String,Any}()
    for name in vSolvedWithInit
        nameAsExpr = Meta.parse(name)
        if haskey(init, nameAsExpr)
            vSolvedWithInitValuesAndUnit[name] = eval( init[nameAsExpr] )
        else
            @warn "Internal issue of TinyModia: $name is assumed to have an init-value, but it is not found."
        end
    end

    # Generate code
    nCrossingFunctions, nAfter, nClocks, nSamples, previousVars, preVars, holdVars = getEventCounters()   
    previousVars = Symbol.(previousVars)
    preVars = Symbol.(preVars)
    holdVars = Symbol.(holdVars)
    
    # Variables added to result
    extraResults = vcat(:time, setdiff([Symbol(u) for u in unknowns], 
                                        Symbol[Symbol(xi_info.x_name_julia)     for xi_info in equationInfo.x_info],
                                        Symbol[Symbol(xi_info.der_x_name_julia) for xi_info in equationInfo.x_info]))
    
    if true # logTiming
#        println("Generate code")
        if useNewCodeGeneration
            @timeit to "generate_getDerivativesNew!" code = generate_getDerivativesNew!(AST, newFunctions, modelModule, equationInfo, [:(_p)], vcat(:time, [Symbol(u) for u in unknowns]), previousVars, preVars, holdVars, :getDerivatives, hasUnits = !unitless)
        else
            @timeit to "generate_getDerivatives!" code = generate_getDerivatives!(AST, equationInfo, [:(_p)], vcat(:time, extraResults), previousVars, preVars, holdVars, :getDerivatives, hasUnits = !unitless)
        end
    else
#        code = generate_getDerivatives!(AST, equationInfo, Symbol.(keys(parameters)), vcat(:time, extraResults), :getDerivatives, hasUnits = !unitless)
        code = generate_getDerivatives!(AST, newFunctions, modelModule, equationInfo, [:(_p)], vcat(:time, extraResults), previousVars, preVars, holdVars, :getDerivatives, hasUnits = !unitless)
    end
    if logCode
        @show mappedParameters
        showCodeWithoutComments(code)
    end

    # Compile code

#    generatedFunction = @RuntimeGeneratedFunction(modelModule, code)
    #getDerivatives = Core.eval(modelModule, code)

    if logTiming
        println("eval code")
        @time @timeit to "eval(code)" Core.eval(modelModule, code)
    else
        Core.eval(modelModule, code)
    end
    getDerivatives = modelModule.getDerivatives

    # If generatedFunction is not packed inside a function, DifferentialEquations.jl crashes
#    getDerivatives(derx,x,m,time) = generatedFunction(derx, x, m, time)

    convertedStartValues = convert(Vector{FloatType}, [ustrip(v) for v in startValues])  # ustrip.(value) does not work for MonteCarloMeasurements
#    @show mappedParameters

#    println("Build SimulationModel")

    model = @timeit to "build SimulationModel" SimulationModel{FloatType, OrderedDict{Symbol,Any}}(modelModule, name, getDerivatives, equationInfo, convertedStartValues, previousVars, preVars, holdVars,
#                                         parameters, vcat(:time, [Symbol(u) for u in unknowns]);
                                         OrderedDict(:(_p) => mappedParameters ), extraResults;
                                         vSolvedWithInitValuesAndUnit, vEliminated, vProperty,
                                         var_name = (v)->string(unknownsWithEliminated[v]),
                                         nz=nCrossingFunctions, nAfter=nAfter,  unitless=unitless)

    # Execute code
    if logExecution
        println("\nExecute getDerivatives")
    #    @show startValues
        derx = deepcopy(convertedStartValues) # To get the same type as for x (deepcopy is needed for MonteCarloMeasurements)
        println("First executions of getDerivatives")
        @timeit to "execute getDerivatives" try
            TimeType = timeType(model)
            @time Base.invokelatest(getDerivatives, derx, convertedStartValues, model, convert(TimeType, 0.0))
            @time Base.invokelatest(getDerivatives, derx, convertedStartValues, model, convert(TimeType, 0.0))
#            @show derx
        catch e
            error("Failed: ", e)
            return nothing
        end
    end
    return model
end

"""
    modelInstance = @instantiateModel(model; FloatType = Float64, aliasReduction=true, unitless=false,
        log=false, logModel=false, logDetails=false, logStateSelection=false, logCode=false, 
        logExecution=logExecution, logCalculations=logCalculations, logTiming=false)

Instantiates a model, i.e. performs structural and symbolic transformations and generates a function for calculation of derivatives suitable for simulation.

* `model`: model (declarations and equations)
* `FloatType`: Variable type for floating point numbers, for example: Float64, Measurements{Float64}, StaticParticles{Float64,100}, Particles{Float64,2000}
* `aliasReduction`: Perform alias elimination and remove singularities
* `unitless`: Remove units (useful while debugging models and needed for MonteCarloMeasurements)
* `log`: Log the different phases of translation
* `logModel`: Log the variables and equations of the model
* `logDetails`: Log internal data during the different phases of translation
* `logStateSelection`: Log details during state selection
* `logCode`: Log the generated code
* `logExecution`: Log the execution of the generated code (useful for timing compilation)
* `logCalculations`: Log the calculations of the generated code (useful for finding unit bugs)
* `logTiming`: Log timing of different phases
* `return modelInstance prepared for simulation`
"""
macro instantiateModel(model, kwargs...)
    modelName = string(model)
    code = :( instantiateModel($model; modelName=$modelName, modelModule=@__MODULE__, source=@__FILE__, $(kwargs...) ) )
    return esc(code)
end

"""
See documentation of macro @instatiateModel
"""
function instantiateModel(model; modelName="", modelModule=nothing, source=nothing, FloatType = Float64, aliasReduction=true, unitless=false,
    log=false, logModel=false, logDetails=false, logStateSelection=false, logCode=false, 
    logExecution=logExecution, logCalculations=logCalculations, logTiming=false, evaluateParameters=false)
    try
        println("\nInstantiating model $modelModule.$modelName")
        resetEventCounters()
        global to = TimerOutput()

        modelStructure = ModelStructure()

        if isexpr(model, :quote)
            model = eval(model) # model defined with macro
        end

        if typeof(model) <: NamedTuple || typeof(model) <: Dict || typeof(model) <: OrderedDict
            if logModel
                @showModel(model)
            end
            if false
                println("drawModel")
                jsonDiagram = drawModel(modelName, model)
                println()
#                println("jsonModel")
#                JSON.print(jsonDiagram, 2)
            end

            if logTiming
                println("Flatten")
                @timeit to "flatten" flattenModelTuple!(model, modelStructure, modelName, to; unitless, log)
            else
                flattenModelTuple!(model, modelStructure, modelName, to; unitless, log)
            end
            pars = OrderedDict{Symbol, Any}([(Symbol(k),v) for (k,v) in modelStructure.parameters])
            if length(modelStructure.equations) > 0
                flatModel = Model() | pars | Model(equations = [removeBlock(e) for e in modelStructure.equations])
            else
                flatModel = Model() | pars 
            end
            printModelStructure(modelStructure, log=false)
            name = modelName
        else
            @show model
            error("Invalid model format")
        end

        allVariables = Incidence[]
        if logTiming
            println("Find incidence")
            @timeit to "findIncidence!" findIncidence!(modelStructure.equations, allVariables)
        else
            findIncidence!(modelStructure.equations, allVariables)
        end
        unique!(allVariables)

        unknowns = setdiff(allVariables, keys(modelStructure.parameters), keys(modelStructure.inputs), [:time, :instantiatedModel, :_leq_mode, :_x])
        Avar, states, derivatives = setAvar(unknowns)
        vActive = [a == 0 for a in Avar]

        if logStatistics
            println("Number of states:    ", length(states))
            if length(unknowns)-length(states) != length(modelStructure.equations)
                println("Number of unknowns:  ", length(unknowns)-length(states))
            end
            println("Number of equations: ", length(modelStructure.equations))
        end
        printArray(modelStructure.parameters, "Parameters:", log=log)
        printArray(modelStructure.inputs, "Inputs:", log=log)
        printArray(modelStructure.outputs, "Outputs:", log=log)
        printArray(states, "Potential states:", log=log)
        printArray(setdiff(unknowns, states), "Unknowns:", log=log)
        printArray(modelStructure.equations, "Equations:", log=log)
        if logDetails
            @show modelStructure.parameters modelStructure.mappedParameters modelStructure.init modelStructure.start modelStructure.variables modelStructure.flows modelStructure.inputs modelStructure.outputs
        end

        unknownsWithEliminated = unknowns
        if aliasReduction
            if log || logTiming
                println("Perform alias elimination and remove singularities")
            end
            if logTiming
                @timeit to "performAliasReduction" equations, unknowns, Avar, G, states, vEliminated, vProperty = performAliasReduction(unknowns, modelStructure.equations, Avar, logDetails, log)
                println("Number of reduced unknowns:  ", length(unknowns)-length(states))
                println("Number of reduced equations: ", length(equations))
            else
                equations, unknowns, Avar, G, states, vEliminated, vProperty = performAliasReduction(unknowns, modelStructure.equations, Avar, logDetails, log)
            end
            printArray(states, "States:", log=log)
            printArray(setdiff(unknowns, states), "Unknowns after alias reduction:", log=log)
            printArray(equations, "Equations after alias reduction:", log=log)
        else
            equations = modelStructure.equations
            G = Vector{Int}[] # Array{Array{Int64,1},1}()
            variablesIndices = OrderedDict(key => k for (k, key) in enumerate(unknowns))

            println("Build incidence matrix")
            unknownsSet = Set(unknowns)
            @timeit to "build incidence matrix" for i in 1:length(equations)
                e = equations[i]
                nameIncidence, coefficients, rest, linear = getCoefficients(e)

                incidence = [] # [variablesIndices[n] for n in nameIncidence if n in unknownsSet]
                for n in nameIncidence
                    if n in keys(variablesIndices)
                        push!(incidence, variablesIndices[n])
                    end
                end
                unique!(incidence)
                push!(G, incidence)
            end
            vEliminated = Int[]
            vProperty   = Int[]
        end

        modStructure = assignAndBLT(equations, unknowns, modelStructure.parameters, Avar, G, states, logDetails, log, logTiming)

        inst = stateSelectionAndCodeGeneration(modStructure, name, modelModule, FloatType, modelStructure.init, modelStructure.start, modelStructure.inputs, modelStructure.outputs,
            vEliminated, vProperty, unknownsWithEliminated, modelStructure.mappedParameters;
            unitless, logStateSelection, logCode, logExecution, logCalculations, logTiming, evaluateParameters)

        if logTiming
            show(to, compact=true)
            println()
        end

        inst #, flatModel

    catch e
        if isa(e, ErrorException)
            println()
            printstyled("Model error: ", bold=true, color=:red)  
            printstyled(e.msg, "\n", bold=true, color=:red) 
            printstyled("Aborting instantiateModel for $modelName in $modelModule\n", bold=true, color=:red)             
            println()
#            Base.rethrow()
        else
            Base.rethrow()
        end
        show(to, compact=true)
        println()
    end
    

end

end
