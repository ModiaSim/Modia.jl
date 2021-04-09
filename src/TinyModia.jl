"""
Main module of TinyModia.

* Developer: Hilding Elmqvist, Mogram AB
* First version: December 2020
* License: MIT (expat)

"""
module TinyModia

using Reexport

export instantiateModel, @instantiateModel

using Base.Meta: isexpr
using DataStructures: OrderedDict

using ModiaBase.Symbolic
using ModiaBase.Simplify
using ModiaBase.BLTandPantelidesUtilities
using ModiaBase.BLTandPantelides
using ModiaBase.Differentiate
using ModiaBase
export ModiaBase

# using RuntimeGeneratedFunctions
# RuntimeGeneratedFunctions.init(@__MODULE__)
@reexport using Unitful
using  Measurements
import MonteCarloMeasurements

include("NamedTupleModels.jl")
include("EvaluateParameters.jl")
include("EventHandler.jl")
include("CodeGeneration.jl")
include("SimulateAndPlot.jl")


# include("IncidencePlot.jl")
# using .IncidencePlot

const drawIncidence = false

const path = dirname(dirname(@__FILE__))   # Absolute path of package directory

const Version = "0.7.3-dev"
const Date = "2021-04-09"

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
        @time performConsistencyCheck(G, Avar, vActive, parameters, unknowns, states, equations, log)
    else
        performConsistencyCheck(G, Avar, vActive, parameters, unknowns, states, equations, log)
    end

    if logTiming
        println("\nIndex reduction (Pantelides)")
        @time assign, Avar, Bequ = pantelides!(G, length(Avar), Avar)
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
            equations[Bequ[i]] = derivative(equations[i], keys(parameters))
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
        @time assign = matching(HG, length(Avar), vActive)
    else
        assign = matching(HG, length(Avar), vActive)
    end
    if logTiming
        println("BLT")
        @time bltComponents = BLT(HG, assign)
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
    Avar = fill(0, length(unknowns))
    states = []
    derivatives = []
    for i in 1:length(unknowns)
        v = unknowns[i]
        if isexpr(v, :call)
            push!(derivatives, v)
            x = v.args[2]
            ix = findfirst(s -> s==x, unknowns)
            Avar[ix] = i
            push!(states, unknowns[ix])
        end
    end
    return Avar, states, derivatives
end

mutable struct ModelStructure
    parameters::OrderedDict{Any,Any}
    mappedParameters::NamedTuple
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

ModelStructure() = ModelStructure(OrderedDict(), (;), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), OrderedDict(), Expr[])

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
        println("  mappedParameters:  ", m.mappedParameters)
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
    merge!(parent.parameters, prependDict(child.parameters, prefix))
    if length(keys(child.parameters)) > 0
#        @show parent.mappedParameters prefix typeof(prefix) child.mappedParameters child.parameters 
        parent.mappedParameters = merge(parent.mappedParameters, (;prefix => (;child.mappedParameters...)) )
#        @show parent.mappedParameters
    end
    merge!(parent.init, prependDict(child.init, prefix))
    if length(keys(child.init)) > 0
#        @show parent.mappedParameters prefix child.init 
        parent.mappedParameters = recursiveMerge(parent.mappedParameters, (;prefix => (;child.mappedParameters...)) )
    end
    merge!(parent.start, prependDict(child.start, prefix))
    if length(keys(child.start)) > 0
        parent.mappedParameters = recursiveMerge(parent.mappedParameters, (;prefix => (;child.mappedParameters...)) )
    end
    merge!(parent.variables, prependDict(child.variables, prefix))
    merge!(parent.flows, prependDict(child.flows, prefix))
    merge!(parent.inputs, prependDict(child.inputs, prefix))
    merge!(parent.outputs, prependDict(child.outputs, prefix))

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

    variablesIndices = OrderedDict(key => k for (k, key) in enumerate(unknowns))

    linearVariables = []
    nonlinearVariables = []
    linearEquations = Vector{Int}()
    nonlinearEquations = []
    G = Vector{Int}[] # Array{Array{Int64,1},1}()
    Gint = Array{Array{Int64,1},1}()
    GcInt = Array{Array{Int64,1},1}()
    for i in 1:length(equations)
        e = equations[i]
        nameIncidence, coefficients, rest, linear = getCoefficients(e)
        incidence = [variablesIndices[n] for n in nameIncidence if n in unknowns]
        unique!(incidence)
        push!(G, incidence)

        if linear && all([c == 1 || c == -1 for c in coefficients]) && rest == 0 && all(n in unknowns for n in nameIncidence)
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

    unique!(linearVariables)
    unique!(nonlinearVariables)
    vSolveInLinearEquations = setdiff(linearVariables, nonlinearVariables)
    if logDetails
        printArray(linearVariables, "linearVariables:", number=false)
        printArray(equations[linearEquations], "linearEquations:", number=false)
#        printArray(nonlinearVariables, "nonlinearVariables:", number=false)
#        printArray(nonlinearEquations, "nonlinearEquations:", number=false)
        @show linearEquations vSolveInLinearEquations Gint GcInt
    end

    GCopy = deepcopy(G)
    AvarCopy = deepcopy(Avar)
    (vEliminated, vProperty, nvArbitrary, redundantEquations) = simplifyLinearIntegerEquations!(GCopy, linearEquations, GcInt, AvarCopy)

    if logDetails
        @show vEliminated vProperty nvArbitrary redundantEquations
    end
    if log
        printSimplifiedLinearIntegerEquations(GCopy, linearEquations, GcInt, vEliminated, vProperty, nvArbitrary, redundantEquations, v->string(unknowns[v]), printTest=false)
    end

    substitutions = OrderedDict()
    if length(vEliminated) - nvArbitrary > 0
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

    replaceLinearIntegerEquations(GCopy[linearEquations], linearEquations, GcInt, unknowns, equations)
#    printArray(equations, "new equations:", log=log)

    reducedEquations = []

    for ei in 1:length(equations)
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
    for i in 1:length(reducedEquations)
        e = reducedEquations[i]
        nameIncidence, coefficients, rest, linear = getCoefficients(e)
        incidence = [reducedVariablesIndices[n] for n in nameIncidence if n in reducedUnknowns]
        unique!(incidence)
        push!(reducedG, incidence)
    end

    reducedAvar, reducedStates, reducedDerivatives = setAvar(reducedUnknowns)
    return reducedEquations, reducedUnknowns, reducedAvar, reducedG, reducedStates, vEliminated, vProperty
end


function stateSelectionAndCodeGeneration(modelStructure, name, modelModule, FloatType, init, start, vEliminated, vProperty, unknownsWithEliminated, mappedParameters;
    unitless=false, logStateSelection=false, logCode=false, logExecution=false, logTiming=false)
    (unknowns, equations, G, Avar, Bequ, assign, blt, parameters) = modelStructure

    function getSolvedEquationAST(e, v)
        (solution, solved) = solveEquation(equations[e], unknowns[v])
        equ = string(equations[e])
        @assert solved "Equation not solved for $(unknowns[v]): $(equations[e])"
#        sol = string(solution)
#        solution = prepend(makeDerivativeVar(solution, components), :m)
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
        solution = makeDerVar(solution, keys(parameters))
        if logExecution
            var = string(unknowns[v])
            return :($solution; println($var, " = ", upreferred.($(solution.args[1]))))
        else
            return solution
        end
    end

    function getResidualEquationAST(e, residualName)
        eq = equations[e] # prepend(makeDerivativeVar(equations[e], components), :m)
        resid = makeDerVar(:(ustrip($(eq.args[2]) - $(eq.args[1]))), keys(parameters))
        residual = :($residualName = $resid)
        if false #logExecution
            return makeDerVar(:(dump($(makeDerVar(eq.args[2]))); dump($(makeDerVar(eq.args[1]))); $residual; println($residualName, " = ", upreferred.(($(eq.args[2]) - $(eq.args[1]))))))
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
#        @show var init[var]
        if var in keys(init)
            value = eval(init[var])
#            value = 0.0
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
#            value = 0.0
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
        @time (success,AST,equationInfo) = getSortedAndSolvedAST(
            G, Array{Array{Int64,1},1}(blt), assign, Avar, Bequ, stateSelectionFunctions;
            unitless, log = logStateSelection, modelName = name)
    else
        (success,AST,equationInfo) = getSortedAndSolvedAST(
            G, Array{Array{Int64,1},1}(blt), assign, Avar, Bequ, stateSelectionFunctions;
            unitless, log = logStateSelection, modelName = name)
    end
#    equationInfo.nz = nCrossingFunctions

    if ! success
        error("Aborting")
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
    if logTiming
        println("Generate code")
#        @time code = generate_getDerivatives!(AST, equationInfo, Symbol.(keys(parameters)), vcat(:time, [Symbol(u) for u in unknowns]), :getDerivatives, hasUnits = !unitless)
        @time code = generate_getDerivatives!(AST, equationInfo, [:(_p)], vcat(:time, [Symbol(u) for u in unknowns]), :getDerivatives, hasUnits = !unitless)
    else
#        code = generate_getDerivatives!(AST, equationInfo, Symbol.(keys(parameters)), vcat(:time, [Symbol(u) for u in unknowns]), :getDerivatives, hasUnits = !unitless)
        code = generate_getDerivatives!(AST, equationInfo, [:(_p)], vcat(:time, [Symbol(u) for u in unknowns]), :getDerivatives, hasUnits = !unitless)
    end
    if logCode
        @show mappedParameters
        showCodeWithoutComments(code)
    end

    # Compile code


#    generatedFunction = @RuntimeGeneratedFunction(modelModule, code)
    #getDerivatives = Core.eval(modelModule, code)
    Core.eval(modelModule, code)
    getDerivatives = modelModule.getDerivatives

    # If generatedFunction is not packed inside a function, DifferentialEquations.jl crashes
#    getDerivatives(derx,x,m,time) = generatedFunction(derx, x, m, time)

    # Execute code
    if logExecution
        println("\nExecute getDerivatives")
        @show startValues
    end
    convertedStartValues = convert(Vector{FloatType}, [ustrip(v) for v in startValues])  # ustrip.(value) does not work for MonteCarloMeasurements
    model = SimulationModel{FloatType}(modelModule, name, getDerivatives, equationInfo, convertedStartValues,
#                                         parameters, vcat(:time, [Symbol(u) for u in unknowns]);
                                         OrderedDict(:(_p) => mappedParameters ), vcat(:time, [Symbol(u) for u in unknowns]);
                                         vSolvedWithInitValuesAndUnit, vEliminated, vProperty,
                                         var_name = (v)->string(unknownsWithEliminated[v]))

    if logExecution
        derx = deepcopy(convertedStartValues) # To get the same type as for x (deepcopy is needed for MonteCarloMeasurements)
#        @time getDerivatives(derx, convertedStartValues, model, convert(FloatType, 0.0))
        Base.invokelatest(getDerivatives, derx, convertedStartValues, model, convert(FloatType, 0.0))

        @show derx
    end
    return model
end

"""
    modelInstance = @instantiateModel(model; FloatType = Float64, aliasReduction=true, unitless=false,
        log=false, logModel=false, logDetails=false, logStateSelection=false, logCode=false, logExecution=false, logTiming=false)

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
* `logExecution`: Log the execution of the generated code (useful for finding unit bugs)
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
    log=false, logModel=false, logDetails=false, logStateSelection=false, logCode=false, logExecution=false, logTiming=false)
    try
        println("\nInstantiating model $modelModule.$modelName")

        modelStructure = ModelStructure()

        if typeof(model) <: NamedTuple
            if logModel
                @showModel(model)
            end
            if logTiming
                println("Flatten")
                @time flattenModelTuple!(model, modelStructure, modelName; unitless, log)
            else
                flattenModelTuple!(model, modelStructure, modelName; unitless, log)
            end
            printModelStructure(modelStructure, log=logDetails)
            name = modelName
        else
            error("Invalid model format")
        end

        allVariables = Incidence[]
        if logTiming
        println("Find incidence")
            @time findIncidence!(modelStructure.equations, allVariables)
        else
            findIncidence!(modelStructure.equations, allVariables)
        end
        unique!(allVariables)
    #    @show allVariables

        unknowns = setdiff(allVariables, keys(modelStructure.parameters), [:time, :instantiatedModel, :_leq_mode])
        Avar, states, derivatives = setAvar(unknowns)
        vActive = [a == 0 for a in Avar]

        if log
            println("Number of states:    ", length(states))
            println("Number of unknowns:  ", length(unknowns)-length(states))
            println("Number of equations: ", length(modelStructure.equations))
        end
        printArray(modelStructure.parameters, "Parameters:", log=log)
        printArray(states, "Potential states:", log=log)
        printArray(setdiff(unknowns, states), "Unknowns:", log=log)
        printArray(modelStructure.equations, "Equations:", log=log)
        if logDetails
            @show modelStructure.parameters modelStructure.mappedParameters modelStructure.init modelStructure.start modelStructure.variables modelStructure.flows
        end

        unknownsWithEliminated = unknowns
        if aliasReduction
            if log || logTiming
                println("Perform alias elimination and remove singularities")
            end
            if logTiming
                @time equations, unknowns, Avar, G, states, vEliminated, vProperty = performAliasReduction(unknowns, modelStructure.equations, Avar, logDetails, log)
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

            for i in 1:length(equations)
                e = equations[i]
                nameIncidence, coefficients, rest, linear = getCoefficients(e)
                incidence = [variablesIndices[n] for n in nameIncidence if n in unknowns]
                unique!(incidence)
                push!(G, incidence)
            end
            vEliminated = Int[]
            vProperty   = Int[]
        end

        modStructure = assignAndBLT(equations, unknowns, modelStructure.parameters, Avar, G, states, logDetails, log, logTiming)

        stateSelectionAndCodeGeneration(modStructure, name, modelModule, FloatType, modelStructure.init, modelStructure.start, vEliminated, vProperty, unknownsWithEliminated, modelStructure.mappedParameters;
            unitless, logStateSelection, logCode, logExecution, logTiming)

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
    end

end

end
