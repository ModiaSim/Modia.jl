"""
Module for structural analysis of models.

* Developer: Hilding Elmqvist, Mogram AB
* First version: July-August 2016
* Copyright (c) 2016-2019: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)

"""
module BasicStructuralTransform

include("../Tearing.jl")

using ..BLTandPantelides
using ..BLTandPantelidesUtilities
using ..Utilities
import ..Instantiation: GetField, This, Der, Symbolic, prettyPrint, prettyfy
import ..Execution: logTiming
import ..Execution: ModiaSimulationModel, extract_results_ida
import ..Execution: Subs, subs
using DataStructures: OrderedDict

using Base.Meta: quot, isexpr
@static if !(VERSION < v"0.7.0-DEV.2005")
    using SparseArrays
end
@static if VERSION < v"0.7.0-DEV.2005"
    notFound = 0
else
    notFound = nothing
end

using DataStructures
using ..Synchronous
using ..SymbolicTransform
# import ..SymbolicTransform: Add

using ..ExactlyRemoveSingularities
using ..StateSelection
import ModiaMath
using ..ModiaLogging

#= Temporarily removed due to problem with PyPlot
using PyPlot
=#
#using ModiaMath.plot

export residue, hide, skew, skewCoords, residue_der
export analyzeStructurally, basicTransformStructurally, setOptions

# For function FDAE
#import Modia
#using MultiBody

const log = false
const logIndexReduction = false
global useIncidenceMatrix = false
global expandArrayIncidence = false
removeSingularitiesDefault = true
global removeSingularities = removeSingularitiesDefault
global tearing = false
global automaticStateSelection = false
const consistencyCheck = true
global newStateSelection = false
global logFDAE = false
global useKinsol = false
global logStatistics

const DebugPrintAssign = true
const indexReduction = true

#=
@show log
@show logIndexReduction
@show useIncidenceMatrix
@show expandArrayIncidence
@show removeSingularities
@show consistencyCheck
@show newStateSelection
@show logFDAE
@show useKinsol
=#

function setOptions(options)
    global removeSingularities = removeSingularitiesDefault
    if haskey(options, :removeSingularities)
        global removeSingularities = options[:removeSingularities]
        @show removeSingularities
        delete!(options, :removeSingularities)
    end

    global tearing = false
    if haskey(options, :tearing)
        global tearing = options[:tearing]
        @show tearing
        delete!(options, :tearing)
    end

    global automaticStateSelection = false
    if haskey(options, :automaticStateSelection)
        global automaticStateSelection = options[:automaticStateSelection]
        @show automaticStateSelection
        delete!(options, :automaticStateSelection)
    end

    global expandArrayIncidence = false
    if haskey(options, :expandArrayIncidence)
        global expandArrayIncidence = options[:expandArrayIncidence]
        @show expandArrayIncidence
        delete!(options, :expandArrayIncidence)
    end

    global useIncidenceMatrix = false
    if haskey(options, :useIncidenceMatrix)
        global useIncidenceMatrix = options[:useIncidenceMatrix]
        @show useIncidenceMatrix
        delete!(options, :useIncidenceMatrix)
    end

    global newStateSelection = false
    if haskey(options, :newStateSelection)
        global newStateSelection = options[:newStateSelection]
        @show newStateSelection
        delete!(options, :newStateSelection)
    end

    global useKinsol = false
    if haskey(options, :useKinsol)
        global useKinsol = options[:useKinsol]
        @show useKinsol
        delete!(options, :useKinsol)
    end

    global logStatistics = true
    if haskey(options, :logStatistics)
        global logStatistics = options[:logStatistics]
        @show logStatistics
        delete!(options, :logStatistics)
    end
end

hide(x) = x

residue(x) = 0 * x
residue_der(x) = 0 * x

skew(x) = [0 -x[3] x[2]; x[3] 0 -x[1]; -x[2] x[1] 0]

skewCoords(R) = [R[3,2], R[1,3], R[2,1]]

# Copied from SymbolicTransform and modified

zero = Float64(0.0)
iZero = Int64(0)
one = Float64(1.0)

function add(e1, e2)
    if e1 == zero || e1 == iZero
        return e2
    elseif e2 == zero || e2 == iZero
        return e1
    else
        Expr(:call, +, e1, e2)
    end
end

function mult(e1, e2)
    if e1 == zero || e2 == zero || e1 == iZero || e2 == iZero
        return zero
    elseif e1 == one
        return e2
    elseif e2 == one
        return e1
    else
        Expr(:call, *, e1, e2)
    end
end

# end copied

findVariables!(variables, states, deriv::Vector, ex) = nothing
# findVariables!(variables, states, deriv::Vector, der::Der) = (push!(states, der.base); push!(deriv, der); nothing)
findVariables!(variables, states, deriv::Vector, s::Symbol) = (
    # if ! (s in [:m, :s, :kg, :N]);
    push!(variables, s);
    if s == :.; push!(variables, s.args[1]); push!(states, s.args[1]); push!(deriv, s.args[2].value) end; nothing)

  function findVariables!(variables, states, deriv::Vector, ex::Expr)
    if isexpr(ex, :ref, 2) && isa(ex.args[1], Number)
        # Skip symbols in unit notations, e.g. 2[kg]
    elseif isexpr(ex, :call, 3) && ex.args[1] == :* && isa(ex.args[2], Number) && isexpr(ex.args[3], :vect)
        # Skip symbols in units on negative numeric literals, eg -2[m/s]
    elseif isexpr(ex, :call)
        if ex.args[1] == :der
            push!(variables, ex.args[2]);
            push!(states, ex.args[2]);
            push!(deriv, ex)
        else
            for arg in ex.args[2:end]  # skip operator name
                findVariables!(variables, states, deriv, arg)
            end
        end
    elseif isexpr(ex, :comparison)
        findVariables!(variables, states, deriv, ex.args[1])
        findVariables!(variables, states, deriv, ex.args[3])
    elseif isexpr(ex, :comprehension)
        locals = []
        omit = [ex.args[2].args[1]]
        findVariables!(locals, states, deriv, ex.args[1])
        findVariables!(locals, states, deriv, ex.args[2])
        for l in locals
            if !(l in omit)
                push!(variables, l)
            end
        end
    else !isexpr(ex, :quote)
        for arg in ex.args
            findVariables!(variables, states, deriv, arg)
        end
    end
    nothing
end

findIncidence!(incidence::Array{Any,1}, ex) = nothing
# findIncidence!(incidence::Array{Any,1}, s::Symbol) = (if ! (s in [:m, :s, :kg, :N]) push!(incidence, s) end; nothing)
findIncidence!(incidence::Array{Any,1}, s::Symbol) = (push!(incidence, s); nothing)
# findIncidence!(incidence::Array{Any,1}, der::Der) = (push!(incidence, der); nothing)
# findIncidence!(incidence::Array{Any,1}, get::GetField) = (push!(incidence, get.name); nothing)
function findIncidence!(incidence::Array{Any,1}, ex::Expr)
    if !isexpr(ex, :quote)
        if ex.head == :call
            if ex.args[1] == hide
            elseif ex.args[1] == Synchronous.previous
                for arg in ex.args[2:end]
                    findIncidence!(incidence, arg)
                end
            elseif ex.args[1] == :der
                push!(incidence, ex)
            else
                for arg in ex.args[2:end]  # skip operator/function name ???
                    findIncidence!(incidence, arg)
                end
            end
        else
            for arg in ex.args
                findIncidence!(incidence, arg)
            end
        end
    end
    nothing
end

function collectIncidence!(equnr, orgEqu, incidenceMatrix, G, DAEvariables, assigned, assignedv)
    # @show assignedv
    for v in G[equnr]
        # @show G[equnr]
        if v == assignedv
        elseif v in DAEvariables
            incidenceMatrix[orgEqu, v] = true
        else
            collectIncidence!(assigned[v], orgEqu, incidenceMatrix, G, DAEvariables, assigned, v)
        end
    end
end

findCoefficients!(coefficients, nonLinear, ex) = nonLinear
findCoefficients!(coefficients, nonLinear, ex::Int64) = (coefficients[:1] = ex; nonLinear)
findCoefficients!(coefficients, nonLinear, ex::Float64) = true
findCoefficients!(coefficients, nonLinear, der::Der) = (coefficients[der] = 1; nonLinear)
findCoefficients!(coefficients, nonLinear, get::GetField) = (coefficients[get.name] = 1; nonLinear)
findCoefficients!(coefficients, nonLinear, s::Symbol) = (coefficients[s] = 1; nonLinear)
function findCoefficients!(coefficients, nonLinear, ex::Expr)
    if !isexpr(ex, :quote) && !isexpr(ex, :line)
        if ex.head in [:(=)]
            for iarg in 1:length(ex.args)
                minus = iarg > 1
                arg = ex.args[iarg]
                coeff = Dict()
                nonLinear = findCoefficients!(coeff, nonLinear, arg) || nonLinear
                for v in keys(coeff)
                    if !minus
                        if !haskey(coefficients, v)
                            coefficients[v] = coeff[v]
                        else
                            coefficients[v] += coeff[v]
                        end
                    else
                        if !haskey(coefficients, v)
                            coefficients[v] = -coeff[v]
                        else
                            coefficients[v] -= coeff[v]
                        end
                    end
                end
            end
        elseif ex.head == :call
            if ex.args[1] in [+, :+]
                for iarg in 2:length(ex.args)
                    arg = ex.args[iarg]
                    coeff = Dict()
                    nonLinear =  findCoefficients!(coeff, nonLinear, arg) || nonLinear

                    if coefficients == nothing || coeff == nothing
                        coefficients = nothing
                    else
                        for v in keys(coeff)
                            if !haskey(coefficients, v)
                                coefficients[v] = coeff[v]
                            else
                                coefficients[v] += coeff[v]
                            end
                        end
                    end
                end
            elseif ex.args[1] in [-, :-]
                for iarg in 2:length(ex.args)
                    minus = length(ex.args) == 2 || iarg > 2
                    arg = ex.args[iarg]
                    coeff = Dict()
                    nonLinear = findCoefficients!(coeff, nonLinear, arg) || nonLinear

                    if coefficients == nothing || coeff == nothing
                        coefficients = nothing
                    else
                        for v in keys(coeff)
                            if !minus
                                if !haskey(coefficients, v)
                                    coefficients[v] = coeff[v]
                                else
                                    coefficients[v] += coeff[v]
                                end
                            else
                                if !haskey(coefficients, v)
                                    coefficients[v] = -coeff[v]
                                else
                                    coefficients[v] -= coeff[v]
                                end
                            end
                        end
                    end
                end
            elseif ex.args[1] in [*, :*] && length(ex.args) == 3
                if typeof(ex.args[2]) == Int64
                    arg = ex.args[3]
                    coeff = Dict()
                    nonLinear = findCoefficients!(coeff, nonLinear, arg) || nonLinear

                    for v in keys(coeff)
                        if !haskey(coefficients, v)
                            coefficients[v] = ex.args[2] * coeff[v]
                        else
                            coefficients[v] += ex.args[2] * coeff[v]
                        end
                    end
                elseif typeof(ex.args[3]) == Int64
                    arg = ex.args[2]
                    coeff = Dict()
                    nonLinear = findCoefficients!(coeff, nonLinear, arg) || nonLinear
                    for v in keys(coeff)
                        if !haskey(coefficients, v)
                            coefficients[v] = ex.args[3] * coeff[v]
                        else
                            coefficients[v] += ex.args[3] * coeff[v]
                        end
                    end
                else
                    nonLinear = true
                end
            else
                nonLinear = true
            end
        elseif ex.head == :block
            for arg in ex.args
                coeff = Dict()
                nonLinear = findCoefficients!(coeff, nonLinear, arg) || nonLinear
            end
        end
    end
    nonLinear
end

function getCoefficients!(coeff, e)
    nonLinear = false
    nonLinear = findCoefficients!(coeff, nonLinear, e)
    for v in keys(coeff)
        if coeff[v] == 0
            delete!(coeff, v)
        end
    end

    if true
        # loglnModia("\nExpression: ", prettyPrint(e))
        if nonLinear
            loglnModia("Equation is nonlinear")

            loglnModia("Coefficients:")
            for v in keys(coeff)
                loglnModia(coeff[v], " * ", prettyfy(v))
            end
        else
            loglnModia("Coefficients:")
            for v in keys(coeff)
                loglnModia(coeff[v], " * ", prettyfy(v))
            end
        end
    end
    nonLinear
end

function testGetCoefficients()
    coeff = Dict()
    nonLinear = getCoefficients!(coeff, :(a + b + a - c - d = 0))
    nonLinear = getCoefficients!(coeff, :(-a + b + a - (c - d) = 0))
    nonLinear = getCoefficients!(coeff, :(-a + b + a - (-5c - (d - e)) = -(-c)))
    nonLinear = getCoefficients!(coeff, :(-a + b + a - 2 * (-5c - (d - e) * 3) = -(-c)))
    nonLinear = getCoefficients!(coeff, :(2 * b - 11b = 0))
    nonLinear = getCoefficients!(coeff, :(2.0 * b - 11b = 0))
    nonLinear = getCoefficients!(coeff, :(a * b = 0))
    nonLinear = getCoefficients!(coeff, :(sin(a) = 0))
    nonLinear = getCoefficients!(coeff, :(a = 5))
    nonLinear = getCoefficients!(coeff, :(-15 + a + 2b = 5))
end


@enum D undefined = 1 F T

mutable struct Status
    statusV::OrderedDict{Int64,D}
    statusG::OrderedDict{Int64,D}
    statusE::OrderedDict{Int64,D}
end

function setVariable!(subsValues, v, val)
    subsValues[GetField(This(), Symbol(v))] = val
    loglnModia("setVariable: ", v, " = ", val)
end

function evalArray(subsValues, vec)
    evalVec = []
    for v in vec
        s = subs(subsValues, v, false)
        push!(evalVec, eval(s))
    end
    evalVec
end

substituteDer(ex) = ex
function substituteDer(ex::Der)
    :((next($(ex.base)) - $(ex.base)) / delta)
end

function substituteDer(ex::Expr)
    if isexpr(ex, :quote)
        ex
    else
        Expr(ex.head, [substituteDer(arg) for arg in ex.args]...)
    end
end

substituteShift(ex) = ex
function substituteShift(ex::GetField)
    :((next($(ex.name))))
end

function substituteShift(ex::Expr)
    if isexpr(ex, :quote)
        ex
    else
        Expr(ex.head, [substituteShift(arg) for arg in ex.args]...)
    end
end

# ------------------------------------------------

function handleSingularities(coefficients, notLinearVariables, unknowns_indices, names, states, statesIndices, nonStateVariables, equations, orgEquIndex, Avar, G, Gsolvable, ESizes, logTiming)
    loglnModia("\nREMOVE SINGULARITIES")
    if log
        @show coefficients notLinearVariables names states
    end
    linearVars = []
    for eCoeff in coefficients
        for v in keys(eCoeff)
            push!(linearVars, v)
        end
    end
    linearVars = unique(linearVars)
    if log
        printSymbolList("\nLinear variables", linearVars, true)
    end
    linearVarsStrings = [string(v) for v in linearVars]

    linearCoefficientMatrix = spzeros(Int64, length(coefficients), length(linearVars))
    for i in 1:length(coefficients)
        eCoeff = coefficients[i]
        for v in keys(eCoeff)
            j = findfirst(isequal(v), linearVars)
            linearCoefficientMatrix[i, j] = eCoeff[v]
        end
    end

    if log
        fullLinearCoefficientMatrix = Array(linearCoefficientMatrix)
        printobj("fullLinearCoefficientMatrix", fullLinearCoefficientMatrix)
    end
    ix = fill(0, 0)
    for i in 1:length(Avar)
        if Avar[i] != 0
            j = findfirst(isequal(names[i]), linearVars)
            if j != notFound
                push!(ix, j)
            end
        end
    end

    notLinearVariables = [i for i in notLinearVariables if i <= length(linearVars)]
    notLinearVariableNames = names[notLinearVariables]
    if log
        @show notLinearVariableNames
    end

    # Build iy, don't include derivatives or potential states
    iy = [i for i in 1:length(linearVars)]
    for v in notLinearVariableNames
        if v in linearVars
            j = findfirst(isequal(v), linearVars)
            iy = setdiff(iy, j)
        end
    end
    iy = setdiff(iy, ix)

    if log
        @show linearCoefficientMatrix ix iy
    end
    if logTiming
        print("Remove singularities:  ")
        @time result = ExactlyRemoveSingularities.removeSingularities(linearCoefficientMatrix, ix, iy)
    else
        result = ExactlyRemoveSingularities.removeSingularities(linearCoefficientMatrix, ix, iy)
    end

    (iya, eqr, ix1, ix2, eqx, A1, A2) = result

    if log
        @show iya eqr ix1 ix2 eqx A1 A2
        printRemoveSingularities(linearCoefficientMatrix, result, linearVarsStrings)
    end

    # Remove redundant equations and constraint equations
    newEquations = []
    newG = Array{Int,1}[]
    newGsolvable = Array{Int,1}[]
    newESizes = []
    eqrOrg = orgEquIndex[eqr]
    eqxOrg = orgEquIndex[eqx]

    equationsUpdated = false
    for i in 1:length(equations)
        e = equations[i]
        g = G[i]
        if !(i in eqrOrg) & !(i in eqxOrg)
            push!(newEquations, e)
            push!(newG, g)
            push!(newGsolvable, Gsolvable[i])
            push!(newESizes, ESizes[i])
        else
            loglnModia("Removed equation: ", prettyPrint(e))
            equationsUpdated = true
        end
    end

    # Add additional equation to set undefined variables
    for i in iya
        eq = :($(GetField(This(), Symbol(linearVarsStrings[i]))) = 0)
        loglnModia("Added equation: ", prettyPrint(eq))
        equationsUpdated = true
        push!(newEquations, eq)
        push!(newG, [findfirst(isequal(linearVars[i]), names)])
        push!(newGsolvable, [findfirst(isequal(linearVars[i]), names)])
        push!(newESizes, size(0))
    end

    # Add new constraint equations
    for i in 1:size(A1, 1)
        e1 = zero
        e2 = zero
        g = []

        for j in 1:size(A1, 2)
            e1 = add(e1, mult(A1[i,j], GetField(This(), Symbol(linearVarsStrings[ix1[j]]))))
            push!(g, findfirst(isequal(linearVars[ix1[j]]), names))
            push!(nonStateVariables, linearVars[ix1[j]])
            printSymbolList("\nNon state variables", nonStateVariables)

            statesIndices = []
            for k in 1:length(states)
                if !(states[k] in nonStateVariables)
                    push!(statesIndices, [unknowns_indices[states[k].name]])
                end
            end
        end

        for j in 1:size(A2, 2)
            e2 = add(e2, mult(A2[i,j], GetField(This(), Symbol(linearVarsStrings[ix2[j]]))))
            push!(g, findfirst(isequal(linearVars[ix2[j]]), names))
        end

        eq = :($e1 + $e2 = 0.0)
        loglnModia("Added equation: ", prettyPrint(eq))
        push!(newEquations, eq)
        push!(newG, g)
        push!(newGsolvable, [])
        push!(newESizes, size(0)) ### Should be generalized
        equationsUpdated = true
    end

    equations = newEquations
    # G::Array{Array{Int,1},1} = newG
    G = newG
    Gsolvable = newGsolvable
    ESizes = newESizes

    if equationsUpdated
        loglnModia("\nUpdated equations")
        for e in equations
            loglnModia(prettyPrint(e))
        end
    end
    return equations, G, Gsolvable, ESizes, nonStateVariables
end


# ---------------------------

#function analyzeStructurally(equations, params, unknowns_indices, deriv, unknownsNames, Avar, statesIndices, states, nonStateVariables, n, realStates, findIncidence!, VSizes, VTypes, ESizes, ETypes, unknowns, partial)

const tSizes = Array{Tuple{Int64,Vararg{Int64,N} where N},1}

function analyzeStructurally(equations, params, unknowns_indices, deriv, unknownsNames, Avar, statesIndices, states, nonStateVariables, n, realStates, findIncidence!, VSizes, VTypes, ESizes, ETypes, unknowns, partial)
  # Build incidence structure as bipartite graph.
    neq = 0
    G = [] # Array{Any}(undef, 0)
    Gsolvable = []
    coefficients = []
    notLinearVariables = []
    orgEquIndex = []

    for eq in equations
        neq += 1
        if DebugPrintAssign loglnModia("\nequation ", neq, ": ", prettyPrint(eq)) end
        incidence = Array{Any,1}()
        findIncidence!(incidence, eq)
        if DebugPrintAssign
            # loglnModia("incidence = ", incidence)
            printSymbolList("incidence", incidence)
        end

        vertices = []
        for inc in incidence
            i = get(unknowns_indices, inc, 0)
            if i == 0
                i = findfirst(isequal(inc), deriv)
                if i != notFound
                    i = i + length(unknownsNames)
                end
            end
            if i != notFound
                push!(vertices, i)
            end
        end

        vertices = unique(vertices)

        if DebugPrintAssign
            loglnModia("vertices = ", vertices)
        end
        push!(G, vertices)

        if removeSingularities || tearing || automaticStateSelection
            # Find linear equations with integer coefficients without offset.
            coeff = Dict()
            solvable = []
            nonLinear = getCoefficients!(coeff, eq)

            @static if VERSION < v"0.7.0-DEV.2005"
                emptySet = []
            else
                emptySet = Set(Any[])
            end
#            @show coeff keys(coeff) keys(params)
#            @show intersect(keys(coeff), keys(params))
            if !nonLinear && (intersect(keys(coeff), keys(params)) == emptySet) && !(1 in keys(coeff))
                push!(coefficients, copy(coeff))
                push!(orgEquIndex, neq)
            else
                for v in vertices
                    push!(notLinearVariables, v)
                end
            end
            for inc in incidence
                if inc in keys(coeff) && coeff[inc] in [1, -1]
                    index = get(unknowns_indices, inc, 0)
                    if index == 0
                        index = findfirst(isequal(inc), deriv)
                        if index != notFound
                            index = index + length(unknownsNames)
                        end
                    end
                    if index != notFound
                        push!(solvable, index)
                    end
                end
            end
            push!(Gsolvable, sort(solvable))
        end
    end

    names = createNames(unknownsNames, Avar)

    if removeSingularities
        (equations, G, Gsolvable, ESizes, nonStateVariables) = handleSingularities(coefficients, notLinearVariables, unknowns_indices, names, states, statesIndices, nonStateVariables, equations, orgEquIndex, Avar, G, Gsolvable, ESizes, logTiming)
    end

    # ------------------------------------------------

    #=
    equationsInfix = []
    for e in equations
        push!(equationsInfix, prettyPrint(e))
    end
    #  equationsInfix = prettyPrint(equations)
    =#
    equationsInfix = string.(prettyPrint.(equations))
    # equationsInfix = makeList(equationsInfix, 1:length(Bequ), Bequ, true)

    # Set sizes of derivative variables
    VSizes = [VSizes; fill("", length(Avar) - length(VSizes))]
    for i in 1:length(Avar)
        if Avar[i] != 0
            VSizes[Avar[i]] = VSizes[i]
        end
    end

    M = length(Avar)

    # ------------------------------------------------

    if consistencyCheck
        loglnModia("\nPANTELIDES CONSISTENCY CHECK")
        printSymbolList("\nUnknowns", names, true, true, Avar)

        loglnModia()

        hequ = []
        for i in 1:length(Avar)
            a = Avar[i]
            if a > 0
                name = names[i]
                der = names[a]
                push!(hequ, "h($name, $der) = 0")
            end
        end
        equationsEG = [equationsInfix; hequ]

        logModia("\nCheck consistency of equations by matching extended equation set: ")
        EG::Array{Array{Int,1},1} = [G; buildExtendedSystem(Avar)]
        EESizes = copy(ESizes)

        #= Double code
        for i in 1:length(Avar)
            a = Avar[i]
            if a > 0
                # equation: v[i] = der(v[i]) => length(equ) = length(v[i])
                push!(EESizes, VSizes[i])
            end
        end
        =#
        for i in length(G) + 1:length(EG)
            push!(EESizes, VSizes[EG[i][1]])
        end

        if partial
            println("Number of additional equations needed: ", length(VSizes) - length(ESizes))
            return nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
        else
            checkSizes(VSizes, EESizes)
        end

        if expandArrayIncidence
            loglnModia("Array handling:")
            vActive = fill(true, length(Avar))
            vLengths = [prod(s) for s in VSizes]
            eLengths = [prod(s) for s in EESizes]

            if log
                @show VSizes
                @show vLengths
                @show eLengths
            end

            if logTiming
                print("Consistency check:     ")
                @time assignEG = matchingArray(EG, vActive, vLengths, eLengths)
            else
                assignEG = matchingArray(EG, vActive, vLengths, eLengths)
            end

            loglnModia("---------------")
            loglnModia("assignEG=")
            for c in assignEG
                loglnModia("  ", c)
            end
        else
            if logTiming
                print("Consistency check:     ")
                @time assignEG = matching(EG, length(Avar))
            else
                assignEG = matching(EG, length(Avar))
            end
        end

        if log
            @show EG
            @show assignEG
        end

        if false
            componentsEG = BLTArray(EG, assignEG, vLengths, eLengths)
            # @show componentsEG
        end

        # checkAssign(assignEG, VSizes, VTypes, EESizes, ETypes, equationsEG, unknownsNames, Avar)

        if !expandArrayIncidence && !all(assignEG .> 0) # Temporary disable check for expandArrayIncidence
            loglnModia("\nThe model is not consistent since all variables could not be assigned:")
            BequEG1 = fill(0, length(EG))
            printUnassigned(equationsEG, unknownsNames, assignEG, Avar, BequEG1)

            loglnModia("\nAssigned equations:")
            printAssignedEquations(equationsEG, unknownsNames, 1:length(EG), assignEG, Avar, fill(0, length(EG)))
            loglnModia()
            error("Translation of model aborted since model not consistent. See log file.")
        else
            loglnModia("Consistency check OK")
        end

        if false # big skip
            if log
                loglnModia("\nAssigned extended equations:")
                printAssignedEquations(equationsEG, unknownsNames, 1:length(EG), assignEG, Avar, fill(0, length(EG)))
            end

            global visiting = []
            #=
            function findUnassigned(G, assigned, e, unassigned, VSizes, ESizes)
                loglnModia("findUnassigned equation $e")
                push!(visiting, e)
                @show visiting
                @show e unassigned
                g = G[e]
                @show g
                int = intersect(unassigned, g )

                if int != []
                    @show int
                    @show e
                    loglnModia("return")
                    return (true, e, int[1])
                else
                    for u in g
                        @show u
                        if ! (u in visiting)
                            (found, eFound, vFound) = findUnassigned(G, assigned, assigned[u], unassigned, VSizes, ESizes)
    #                       pop!(visiting)
                            if found
                                @show found
                                return (true, eFound, vFound)
                            end
                        end
                    end
                end
            end
            =#

        end # big skip
    end

    Bequ = fill(0, length(G))

    # ------------------------------------------------

    # Add function arguments and unnest function.
    function reduceDAEIndex()
        loglnModia("\nINDEX REDUCTION")
        if expandArrayIncidence
            vLengths = [prod(s) for s in VSizes]
            eLengths = [prod(s) for s in ESizes]

            if log
                @show VSizes
                @show vLengths
                @show eLengths
            end

            if logTiming
                print("PantelidesArray:         ")
                @time (assign, Avar, Bequ) = pantelidesArray!(G, length(Avar), Avar, vLengths, eLengths)
            else
                (assign, Avar, Bequ) = pantelidesArray!(G, length(Avar), Avar, vLengths, eLengths)
            end

            if log
                @show Avar Bequ
            end

            loglnModia("---------------")
            loglnModia("assign=")

            for c in assign
                loglnModia("  ", c)
            end

        else
            if logTiming
                print("Pantelides:            ")
                @time (assign, Avar, Bequ) = pantelides!(G, length(Avar), Avar)
            else
                (assign, Avar, Bequ) = pantelides!(G, length(Avar), Avar)
            end
        end

        # checkSizes(VSizes, ESizes)   # VSizes, ESizes not up-to-date

        if logIndexReduction
            # @show assign
            # @show Avar
            # @show Bequ
            if !expandArrayIncidence
                loglnModia("\nAssigned equations after index reduction:")
                printAssignedEquations(equationsInfix, unknownsNames, 1:length(G), assign, Avar, Bequ)
            end
        end

        # Commented in order to be able to expandArrayIncidence
        # printUnassigned(equationsInfix, unknownsNames, assign, Avar, Bequ)
        loglnModia()
        ###  checkAssign(assign, VSizes, VTypes, ESizes, ETypes, equationsInfix, unknownsNames, Avar)

        #=
        components = BLT(G, assign)

        if log
            # @show components

            loglnModia("\nSorted equations:")
            printSortedEquations(equationsInfix, unknownsNames, components, assign, Avar, Bequ)
        end
        =#
        #=
        AG = [G; buildFullIncidence(length(Avar)-length(Bequ), length(Avar))]
        assignAG = matching(AG, length(Avar))
        componentsAG = BLT(AG, assignAG)
        if false # log
            loglnModia("\nAugmented system.")
            @show AG
            @show assignAG
            @show componentsAG

            loglnModia("\nSorted augmented equations:")
            equationsAG = [equationsInfix; fill("der(...)", length(Bequ)-length(equations)); fill("full", length(Avar)-length(Bequ))]
            BAG = [Bequ; fill(0, length(Avar)-length(Bequ))] # Additional equations not differentiated
            printSortedEquations(equationsAG, unknownsNames, componentsAG, assignAG, Avar, BAG)

            if ! all(assignAG .> 0)
                loglnModia("The model is not consistent since all variables could not be assigned:")
                printUnassigned(equationsAG, unknownsNames, assignAG, Avar, BAG)
                loglnModia()
            else
                loglnModia("Augmented system OK")
            end
        end
        =#
    end # reduceDAEIndex

    # ----------------------------

    assign = []
    if indexReduction
        reduceDAEIndex()
    end

    loglnModia("\nSTATE SELECTION")
    IG = Array{Array{Int64,1},1}(G)

    for i in 1:length(IG)-length(Gsolvable)
        push!(Gsolvable, [])  # This is conservative since highest derivative is not marked as solvable.
    end

    if automaticStateSelection
        startValues = []
        fixedFlags = []
        for (name, var) in unknowns
            push!(startValues, var.start)
            push!(fixedFlags, var.fixed)
        end
        @show startValues
        @show fixedFlags
        components = Array{Array{Int64,1},1}(BLT(IG, assign))
        vNames = makeList(unknownsNames, 1:length(assign), Avar) # ::Vector{String}
        eqGraph = StateSelection.getSortedEquationGraph(IG, Gsolvable, components, assign, Avar, Bequ, vNames, withStabilization=false)
    end

    if newStateSelection
        # Reduce graph
        newG = copy(IG)
        newAssign = copy(assign)
        # newG = [newG[i] for i in 1:length(newG) if Bequ[i] == 0]
    end

    vActive = fill(true, length(Avar))
    if ! automaticStateSelection
        vActive[statesIndices] .= false
    else
        stateIndices = eqGraph.Vx[ findall(x->x>0, eqGraph.Vx) ]   # Only use Vx indices that are > 0 (= no lambda variables)
        @show names Avar stateIndices
        @show Gsolvable
        for i in 1:length(stateIndices)
            j = stateIndices[i]
            if j in Avar  # selected state is a derivative, but not a lambda variable and not an algebraic variable
                # search for corresponding variable
                alias = 0
                for ie in 1:length(G)
                    g = Gsolvable[ie]
                    if j in g && length(g) == 2 && sort(g) == sort(G[ie])
                        alias = if g[1] == j; g[2] else g[1] end
                        break;
                    end
                end
                @show alias i
                stateIndices[i] = alias
            end
        end

        @show stateIndices
        vActive[stateIndices] .= false
        stateNames = [replace(replace(string(name), "(" => "_"), ")" => "") for name in names[stateIndices]]
        @show stateNames
        @show realStates
        realStates = [GetField(This(), Symbol(name)) for name in stateNames]
        setRealStates(realStates)
    end

    if log
        println("\nNot active:")
        for i in 1:length(vActive)
            if !vActive[i]
                @show i
            end
        end
        @show vActive
    end
    if expandArrayIncidence
        # @show IG
        if logTiming
            print("MatchingArray:           ")
            @time assignIG = matchingArray(IG, vActive, vLengths, eLengths)
        else
            assignIG = matchingArray(IG, vActive, vLengths, eLengths)
        end

        if log
            @show assignIG
        end
        loglnModia("---------------")
        loglnModia("assignIG=")

        for c in assignIG
            loglnModia("  ", c)
        end
        # assignIGTemp = [[e[1] for e in c][1] for c in assignIG]  # Not perfect within strong components

    else
        if logTiming
            print("Matching:              ")
            @time assignIG = matching(IG, length(Avar), vActive)
        else
            assignIG = matching(IG, length(Avar), vActive)
        end
        assignIGTemp = assignIG
    end
    ###  checkAssign(assignIG, VSizes, VTypes, ESizes, ETypes, equationsInfix, unknownsNames, Avar, vActive)

    # @show assignIG
    if !expandArrayIncidence
        (invAssign, unAssignedVariables) = invertAssign(assignIGTemp, length(Bequ))  ### Enabled again
    end

    unassignedNames = [] # unknownsNames[unAssignedVariables]

    # Check number of variables
    variableEquationDifference = length(IG) + length(realStates) - length(assignIG)
    if variableEquationDifference != 0
       error("\nThe number of equations and the number of unknowns/states does not match.",
             "\n   Number of equations: ", length(IG),
             "\n   Number of variables: ", length(assignIG),
             "\n   Number of continuous states: ", length(realStates),
             "\nNote: Number of equations - Number of variables + Number of continuous states must be zero.",
             variableEquationDifference > 0 ?
                 "\nIt might be that " * string(variableEquationDifference) * " states must be defined to be non-states" *
                 "\n(note, Modia does not yet support automatic state selection)." : "")
    end

    if logStatistics
        println("Number of equations: ", length(IG))
        println("Number of variables: ", length(assignIG))
        println("Number of continuous states: ", length(realStates))

        if length(nonStateVariables) > 0
            println("Number of non states: ", length(nonStateVariables))
        end

        # printSymbolList("Unknowns: ", unknownsNames)
        # printSymbolList("Real states: ", realStates)
        # printSymbolList("Non states: ", nonStateVariables)
    end

    if false # log
        @show IG
        @show assignIG
        @show invAssign
        @show Avar
        @show Bequ

        loglnModia("\nAssigned equations after index reduction and state selection:")
        printAssignedEquations(equationsInfix, unknownsNames, 1:length(G), assignIGTemp, Avar, Bequ)

        printUnassigned(equationsInfix, unknownsNames, assignIGTemp, Avar, Bequ, vActive)
        loglnModia()
    end

    loglnModia("\nBLT TRANSFORMATION")
    if expandArrayIncidence
        loglnModia("IG:")
        vNames = makeList(unknownsNames, 1:length(assignIG), Avar) # ::Vector{String}
        ii = 0
        for e in IG
            loglnModia("  ", e)
            ii += 1
            logModia(ii, " ")
            for v in e
                logModia(vNames[v], " ")
            end
            loglnModia()
        end
        # @show vNames[statesIndices]

        # @show IG
        # @show assignIG vLengths eLengths             ")
        if logTiming
            print("BLTArray:                ")
            @time componentsIG = BLTArray(IG, assignIG, vLengths, eLengths)
        else
            componentsIG = BLTArray(IG, assignIG, vLengths, eLengths)
        end

        # (invAssign, unAssignedVariables) = invertAssign(assignIG, length(Bequ))
        equationsIG = [equationsInfix; fill("der(...)", length(Bequ) - length(equations))]
        for e in 1:length(Bequ)
            if Bequ[e] > 0
                equationsIG[Bequ[e]] = "DER(" * string(prettyPrint(equationsIG[e])) * ")"
            end
        end
        (invAssign, assignArray) = compactify(componentsIG, assignIG, vLengths, eLengths, vNames, equationsIG)
        # @show invAssign

        loglnModia()
        # @show componentsIG
        loglnModia("componentsIG=")
        for c in componentsIG
            loglnModia("  ", c)
        end
        loglnModia("Compacting")
        # All array variable elements belong to the same strongly connected components. So it's sufficient to pick out just the first element of the tuple. And remove duplicates.
        loglnModia()
        componentsIGorg = componentsIG
        # @show componentsIGorg
        componentsIG = [unique([e[1] for e in c]) for c in componentsIG]
        if log
            @show componentsIG
        end
        loglnModia("---------------")
        loglnModia("assignIG=")
        for c in assignIG
            loglnModia("  ", c)
        end

        if log
            @show assignIG
        end
        assignIGorg = assignIG
        # @show assignIGorg
        println()
        # assignIG = [[e[1] for e in c][1] for c in assignIG]  # Not perfect within strong components
        assignIG = assignArray
        println()
        # @show assignIG
        loglnModia("CHECK:")

        if log
            @show assignIG
        end
    else
        if logTiming
            print("BLT:                   ")
            @time componentsIG = BLT(IG, assignIG)
        else
            componentsIG = BLT(IG, assignIG)
        end

        if tearing
            loglnModia("\nTearing:")

            equationsIG = [equationsInfix; fill("der(...)", length(Bequ) - length(equations))]

            if true # log
                loglnModia("\nSorted state equations before tearing:")
                printSortedEquations(equationsIG, unknownsNames, componentsIG, assignIG, Avar, Bequ)
                loglnModia()
            end

            tornComponents = []
            for c in componentsIG
                if length(c) == 1
                    push!(tornComponents, c)
                else
                    es = Array{Int64,1}(c)
                    vs = invAssign[es]

                    td = TraverseDAG(IG, length(assignIG))
                    (eSolved, vSolved, eResidue, vTear) = tearEquations!(td, Gsolvable, es, vs)

                    for e in eSolved
                        push!(tornComponents, [e])
                    end

                    push!(tornComponents, eResidue)

                    for i in 1:length(eSolved)
                        assignIG[vSolved[i]] = eSolved[i]
                    end
                    for i in 1:length(eResidue)
                        # Mark by 0 that the equation can not be solved symbolically even if possible.
                        assignIG[vTear[i]] = 0 # eResidue[i]
                    end

#                   (invAssign, unAssignedVariables) = invertAssign(assignIG, length(Bequ))

                    if length(c) > length(eResidue)
                        loglnModia("Reduced system of equation size from $(length(c)) to $(length(eResidue))")
                    else
                        loglnModia("No reduction of system of equation size $(length(c))")
                    end
                end
            end
            componentsIG = tornComponents
        end
    end

    if log
        @show componentsIG
    end

    if newStateSelection
        (sortedEquations, equation_code) = performNewStateselection(unknownsNames, deriv, equations, equationsInfix, newG, assignIG, newAssign, Avar, Bequ, expandArrayIncidence)
    else
        sortedEquations = nothing
        equation_code = nothing
    end

    if newStateSelection || useKinsol # && false ## Disable simulation for 0.7
        generateCode(newStateSelection, useKinsol, params, realStates, unknownsNames, deriv, equations, componentsIG, assignIG, Avar, Bequ, VSizes, ESizes, invAssign, sortedEquations, equation_code, unknowns)

        return nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
    end

    # ---
    incidenceMatrix = spzeros(Bool, length(IG), length(assignIG))

    for i in 1:length(IG)
        for j in 1:length(IG[i])
            incidenceMatrix[i, IG[i][j]] = true
        end
    end

    if false # PrintFinalIncidenceMatrix
        intIncidenceMatrix = [if incidenceMatrix[i,j]; 1 else 0 end for i in 1:size(incidenceMatrix, 1), j in 1:size(incidenceMatrix, 2)]
        @show intIncidenceMatrix
    end

    DAEvariables = []
    for i in 1:n
        if !vActive[i] || i > length(unknownsNames)  # state or derivative
            push!(DAEvariables, i)
        end
    end

    # @show DAEvariables

    derVar = length(unknownsNames) + 1:n
    # @show derVar
    derEqu = assignIG[derVar]
    # @show derEqu

    if useIncidenceMatrix
        incidenceMatrix = spzeros(Bool, length(IG), length(assignIG))

        for equnr in derEqu
            collectIncidence!(equnr, equnr, incidenceMatrix, IG, DAEvariables, assignIG, 0)
        end
        # incidenceMatrix = incidenceMatrix[derEqu, DAEvariables]
        # @show incidenceMatrix

        n2 = div(length(DAEvariables), 2)
        # Treat state and derivative as same incidence:
        incidenceMatrix = incidenceMatrix[derEqu, DAEvariables[1:n2]] .| incidenceMatrix[derEqu, DAEvariables[n2 + 1:length(DAEvariables)]]

        if false # PrintFinalIncidenceMatrix
            intIncidenceMatrix = [if incidenceMatrix[i,j]; 1 else 0 end for i in 1:size(incidenceMatrix, 1), j in 1:size(incidenceMatrix, 2)]
            @show intIncidenceMatrix
        end

        intIncidenceMatrix = nothing
    end

    equationsIG = [equationsInfix; fill("der(...)", length(Bequ) - length(equations))]

    if true # log
        loglnModia("\nSorted state equations:")
        printSortedEquations(equationsIG, unknownsNames, componentsIG, assignIG, Avar, Bequ)
        loglnModia()
    end

    variables = [unknownsNames; deriv]
    # solveSortedEquations(equationsIG, [unknownsNames; deriv], componentsIG, assignIG, Avar, Bequ)

    if false # ! all([assignIG[i] > 0 || ! vActive[i] for i in 1:length(vActive)]) # Temporary disable check for expandArrayIncidence
        loglnModia("\nThe model is not consistent since all variables could not be assigned:")

        @show assignIG
        @show Avar
        @show Bequ
        error("Translation of model aborted")
    end

    loglnModia("End of structural processing.")

    return equations, equationsIG, variables, assignIG, componentsIG, Avar, Bequ, states, deriv, unassignedNames, incidenceMatrix, VSizes, VTypes, ESizes, ETypes
end

# -----------------------------------------------------------------------------

"""
Experimental code for new state selection approach
Add function arguments
"""
function  performNewStateselection(unknownsNames, deriv, equations, equationsInfix, newG, assignIG, newAssign, Avar, Bequ, expandArrayIncidence)
    loglnModia("\nNEW STATE SELECTION")

    # Remove non-leading variables from IG
    if log
        @show IG Avar
    end
    (orgBIndex, derOrder) = invertDer(Bequ)
    # @show orgBIndex
    # IG = [[j for j in G[i] if Bequ[i] !=0 || Avar[j] == 0] for i in 1:length(IG)]
    # IG = [[j for j in G[i] if Bequ[i] != 0 || orgBIndex[i] != 0 && !(j in IG[orgBIndex[i]])] for i in 1:length(IG)]
    # Gsolvable = fill([], length(IG))  # Temporary
    # Gsolvable = newG

    Gsolvable = []
    vNames = makeList(unknownsNames, 1:length(assignIG), Avar) # ::Vector{String}
    @show vNames

    for eq in equations
        coeff = Dict()
        nonLinear = false
        nonLinear = findCoefficients!(coeff, nonLinear, eq)
        solvable = []

        for (var, value) in coeff
            if value in [1, -1]
                ind = findfirst(isequal(string(var)), vNames)
                if ind != notFound
                    push!(solvable, ind)
                end
            end
        end
        push!(Gsolvable, sort(solvable))
    end

    for i in length(Gsolvable) + 1:length(newG)
        solvable = [Avar[v] for v in newG[i] if Avar[v] > 0]
        push!(Gsolvable, solvable)
    end

    if log
        @show length(newAssign) length(componentsIG)
        @show length(newG) newG length(Avar) length(Bequ)
    end

    if expandArrayIncidence
        if log
            @show IG
            @show newAssign vLengths eLengths
        end

        print("BLTArray:                ")
        @time bltA = BLTArray(newG, newAssign, vLengths, eLengths)
        @show bltA
        blt = [unique([e[1] for e in c]) for c in bltA]
        @show blt
    else
        print("BLT     :                ")
        @time blt = BLT(newG, newAssign)
    end

    if log
        @show componentsIG
    end

    blt = [Array{Int64,1}(c) for c in blt]

    printSymbolList("\nvNames", vNames, true, true, Avar)
    loglnModia()

    if log
        @show newG Gsolvable blt newAssign Avar Bequ
    end

    newAssign = [[e[1] for e in c][1] for c in newAssign]  # Not perfect within strong components

    if log
        @show newAssign
        println("\ngetSortedEquationGraph:")
        @show newG Gsolvable blt newAssign Avar Bequ vNames
        println("sortedEquations = getSortedEquationGraph(newG, Gsolvable, blt, newAssign, Avar, Bequ, vNames)")
    end

    sortedEquations = getSortedEquationGraph(newG, Gsolvable, blt, newAssign, Avar, Bequ, vNames)
    loglnModia()

    if log
        @show sortedEquations.isorted
        @show sortedEquations.irbeg
        @show sortedEquations.ESorted
        @show sortedEquations.ESolved
        @show sortedEquations.rcat
        @show sortedEquations.Vx
        @show sortedEquations.VxRev
        @show sortedEquations.Vderx
        @show sortedEquations.VderxRev
        @show sortedEquations.nc
        @show sortedEquations.Er
        @show sortedEquations.ider0n2
        @show sortedEquations.ider1n1
    end

    printSortedEquationGraph(sortedEquations)

    loglnModia("\nSorted variables")
    for i in sortedEquations.ESolved
        if i > 0
            loglnModia(i, ": ", vNames[i])
        else
            loglnModia("  : r[", string(-i), "]")
        end
    end

    loglnModia("\nSorted equations")
    printAssignedEquations(equationsInfix, unknownsNames, sortedEquations.ESorted, fill(0, length(vNames)), Avar, Bequ)

    variables = [unknownsNames; deriv]
    equation_code = newDifferentiateEquations(equations, variables, Avar, Bequ, sortedEquations.ESorted, sortedEquations.ESolved)
    return sortedEquations, equation_code
end

# -----------------------------------------------------------------------------

"""
Experimental code for code generation with Kinsol. Should be merged.
"""
function generateCode(newStateSelection, useKinsol, params, realStates, unknownsNames, deriv, equations, componentsIG, assignIG, Avar, Bequ, VSizes, ESizes, invAssign, sortedEquations, equation_code, unknowns)
    func = "\nmodule ModuleFDAE\n"
    func *= "import Modia\n"
    # func *= "import MultiBody\n"
    func *= "import ModiaMath\n"
    func *= "\nfunction FDAE(_mod, t0, _x, _der_x, _res, _w)\n"

    func *= "\n"
    func *= "  # Function to solve systems of equations\n"

    func *= "  function solve(y, f!)
    info = ModiaMath.NonlinearEquationsInfo(\"\", length(y), f!)
    ModiaMath.solveNonlinearEquations!(info, y)
    return y
  end" * "\n\n"

    if logFDAE
        func *= "  @show _x _der_x\n"
    end

    func *= "  # Set parameters\n"

    for (p, v) in params
        pName = replace(string(p), "." => "_")
        func *= "  " * pName * " = " * string(v) * "\n"
        if logFDAE
            func *= "  " * "@show $pName" * "\n"
        end
    end
    func *= "\n"

    if newStateSelection
        vNames = makeList(unknownsNames, 1:length(assignIG), Avar)

        for i in 1:length(sortedEquations.Vx)
            xName = replace(string(vNames[sortedEquations.Vx[i]]), "." => "_")
            func *= "  $(xName) = _x[$i]\n"
            if logFDAE
                func *= "  " * "@show $xName" * "\n"
            end
        end

        for i in 1:length(sortedEquations.Vderx)
            if sortedEquations.Vderx[i] > 0
                der_xName = replace(string(vNames[sortedEquations.Vderx[i]]), "." => "_")
                func *= "  $(der_xName) = _der_x[$i]\n"
                if logFDAE
                    func *= "  " * "@show $der_xName" * "\n"
                end
            end
        end

        for i in eachindex(sortedEquations.ider0n2)
            iderx = sortedEquations.VxRev[ sortedEquations.ider0n2[i] ]
            ix    = sortedEquations.VxRev[ sortedEquations.ider1n1[i] ]
            func *= "  _r[$i]   = _der_x[$iderx] - _x[$ix]\n"
        end
        # func *= "  _r[1] = _der_x[1] - _x[2]\n"
        # func *= "  " * "println(\"next\")" * "\n"
        for (solution, solved) in equation_code
            pp = prettyPrint(solution)
            pp = replace(string(pp), "der(" => "(der_")
            pp = replace(string(pp), "der_der_" => "der2_")
            if logFDAE
                func *= "  " * "println(\"$pp\")" * "\n"
            end
            func *= "  " * pp * "\n"
        end
    else # Textual FDAE
#      println("\nGENERATING DAE FUNCTION")
        if false # log
            loglnModia("\nSorted state equations:")
            printSortedEquations(equations, unknownsNames, componentsIG, assignIG, Avar, Bequ)
            loglnModia()
        end

        func *= "  time = t0\n"
        
        func *= "  # Copy state vector and derivative vector to each state and derivative\n"
        for i in 1:length(realStates)
            xName = replace(string(realStates[i]), "." => "_")
            xName = replace(xName, "this_" => "")
            derName = "der_" * xName
            func *= "  $(xName) = _x[$i]\n"
            func *= "  $(derName) = _der_x[$i]\n"

            if logFDAE
                func *= "  " * "@show $xName" * "\n"
                func *= "  " * "@show $derName" * "\n"
            end
        end

        variables = [unknownsNames; deriv]
        solvedComponents, algebraic = differentiateSortedEquations(equations, variables, componentsIG, assignIG, Avar, Bequ)
        vNames = makeList(unknownsNames, 1:length(assignIG), Avar)

        # Set sizes of derivative variables ### Duplication of code
        VSizes = [VSizes; fill("", length(Avar) - length(VSizes))]
        for i in 1:length(Avar)
            if Avar[i] != 0
                VSizes[Avar[i]] = VSizes[i]
            end
        end

        ESizes = [ESizes; fill("", length(Bequ) - length(ESizes))]
        for i in 1:length(Bequ)
            if Bequ[i] != 0
                ESizes[Bequ[i]] = ESizes[i]
            end
        end

        #=       func *= "    ("
        for v in vNames
            func *= v
            func *= ", "
        end
        func *= ") = _x" * "\n"
        =#

        # Traverse all BLT blocks
        func *= "\n"
        func *= "  # Equations\n"
        globalResidualIndex = 0
        for blockIndex in 1:length(solvedComponents)
            c = solvedComponents[blockIndex]
            residualIndex = 0
            if length(c) > 1
                func *= "\n"
                func *= "  # Code for solving systems of equations\n"
                # Generate code for solving systems of equations
                # m = length(c)
                # Calculate sizes of y vecor
                m = 0
                for j in 1:length(componentsIG[blockIndex])
                    equIndex = componentsIG[blockIndex][j]
                    varIndex = invAssign[equIndex]
                    if varIndex > 0    ### Should not happen
                        m += prod(VSizes[varIndex])
                    else
                        println("                 NOT ASSIGNED:")
                        @show equIndex
                        println(equationsInfix[equIndex])
                        println()
                    end
                end

                func *= "  _y = solve(ones($m), "
                func *= " function (info, _y, _r)" * "\n"
                if logFDAE
                    func *= "    println(\"f!\")\n"
                end
                # (assignedVar, unAssigned) = invertAssign(assignIG)
                assignedVar = invAssign

                vIndex = 0
                func *= "    # Copy unknowns vector to each unknown\n"
                for j in 1:length(componentsIG[blockIndex])
                    equIndex = componentsIG[blockIndex][j]
                    # @show equIndex
                    varIndex = assignedVar[equIndex]
                    # @show varIndex

                    #=
                    (orgIndexVar, derOrderVar) = invertDer(Avar)
                    if varIndex == 0
                        @show blockIndex
                        @show componentsIGorg[blockIndex] componentsIG[blockIndex]
                        @show j
                        @show equIndex
                        @show assignedVar
                        @show assignedVar[equIndex]
                        @show orgIndexVar
                    end
                    vName = variables[orgIndexVar[varIndex]] ### [1]
                    vName = replace(string(vName), "." => "_")
                    vName = replace(vName, "this_" => "")
                    if derOrderVar[varIndex] == 1  ### Fix for higher order needed
                        vName = "der_"*vName
                    elseif derOrderVar[varIndex] == 2  ### Fix for higher order needed
                        vName = "der2_"*vName
                    end
                    =#
                    if varIndex > 0

                        vName = vNames[varIndex]
                        vName = replace(string(vName), "." => "_")

                        if VSizes[varIndex] == ()
                            vIndex += 1
                            func *= "    $(vName) = _y[$vIndex]" * "\n"
                        else
                            i1 = vIndex + 1
                            s = prod(VSizes[varIndex])
                            i2 = vIndex + s
                            func *= "    $(vName) = reshape(_y[$i1:$i2], $(VSizes[varIndex]))" * "\n"
                            vIndex += s
                            # println(func)
                        end

                        if logFDAE
                            func *= "    " * "@show $vName" * "\n"
                        end
                    end
                end
                func *= "\n"
                func *= "    # Create residuals for system of equations\n"
            end

            # for (e, solved) in c
            for j in 1:length(componentsIG[blockIndex])
                (e, solved) = c[j]
                equIndex = componentsIG[blockIndex][j]
                varIndex = invAssign[equIndex]

                if !solved && length(c) > 1
                    if ESizes[equIndex] == ()
                        residualIndex += 1
                        e = :(_r[$residualIndex] = $(e.args[1]) - $(e.args[2]))
                    else
                        i1 = residualIndex + 1
                        s = prod(ESizes[equIndex])
                        i2 = residualIndex + s
                        e = :(_r[$i1:$i2] = reshape($(e.args[1]) - $(e.args[2]), ($s,)))
                        residualIndex += s
                    end
                elseif !solved
                    println("useKinsol")
                    println(e)
                    # globalResidualIndex += 1
                    # e = :(_res[$globalResidualIndex] = $(e.args[1]) - $(e.args[2]))
                end
                pp = prettyPrint(e)
                pp = replace(string(pp), "der(" => "(der_")
                pp = replace(string(pp), "der_der_" => "der2_")

                if logFDAE
                    if length(c) > 1
                        func *= "  "
                    end
                    func *= "  " * "println(\"$pp\")" * "\n"
                end

                if length(c) > 1
                    func *= "  "
                end
                func *= "  " * string(pp) * "\n"
            end

            if length(c) > 1
                # m = length(c)
                func *= "  end)" * "\n"
                # func *= "  _y = solve(f!, ones($m))" * "\n"
                if logFDAE
                    func *= "  @show _y" * "\n"
                end
                func *= "\n"

                func *= "  # Copy unknowns vector to each unknown\n"
                vIndex = 0
                for j in 1:length(componentsIG[blockIndex])
                    equIndex = componentsIG[blockIndex][j]
                    varIndex = assignedVar[equIndex]
                    #=
                    (orgIndexVar, derOrderVar) = invertDer(Avar)
                    vName = variables[orgIndexVar[varIndex]] ### [1]
                    vName = replace(string(vName), "." => "_")
                    vName = replace(vName, "this_" => "")
                    if derOrderVar[varIndex] == 1  ### Fix for higher order needed
                        vName = "der_"*vName
                    elseif derOrderVar[varIndex] == 2  ### Fix for higher order needed
                        vName = "der2_"*vName
                    end
                    =#
                    if varIndex > 0

                        vName = vNames[varIndex]
                        vName = replace(string(vName), "." => "_")

                        if VSizes[varIndex] == ()
                            vIndex += 1
                            func *= "    $(vName) = _y[$vIndex]" * "\n"
                        else
                            i1 = vIndex + 1
                            s = prod(VSizes[varIndex])
                            i2 = vIndex + s
                            func *= "    $(vName) = reshape(_y[$i1:$i2], $(VSizes[varIndex]))" * "\n"
                            vIndex += s
                            # println(func)
                        end

                        if logFDAE
                            func *= "    " * "@show $vName" * "\n"
                        end
                    end
                end
                func *= "\n"
            end
        end

        func *= "\n"
        func *= "  # Create residual vector for derivatives\n"

        for i in 1:length(realStates)
            xName = replace(string(realStates[i]), "." => "_")
            xName = replace(xName, "this_" => "")
            derName = "der_" * xName
            globalResidualIndex += 1
            func *= "  " * "_res[$globalResidualIndex] = _der_x[$i] - $derName" * "\n"
        end

        n = length(realStates)
        x0 = zeros(n)
        der_x0 = zeros(n)
        r = zeros(n)
        nc = n
    end

    if logFDAE
        func *= "  @show _res" * "\n"
    end
    func *= "end\n"
    func *= "\nend\n"

    loglnModia(func)

    # Allow editing of FDAE
    if true
        open("FDAE.jl", "w") do f
        write(f, func)
      end
    end

    if false
        open("FDAE.jl") do f
        func = readstring(f)
      end
    end

    FUNC = Meta.parse(func)
    fdae = eval(FUNC)

    if newStateSelection
        loglnModia("Initial conditions")
        x0 = fill(0.0, length(sortedEquations.Vx))
        @show vNames[sortedEquations.Vx]
        for (name, var) in unknowns
            if var.start != nothing
                loglnModia(name, " = ", var.start)
                k = findfirst(isequal(string(name)), vNames[sortedEquations.Vx])
                if k != notFound
                    x0[k] = var.start  # Need to take care of vectors.
                else
                    println("Not found: $(string(name)), start=$(var.start)")
                end
            end
        end

        @show x0
        der_x0 = zeros(length(sortedEquations.Vderx))
        r = zeros(length(sortedEquations.Er))
        @show r
        nc = sortedEquations.nc
    end

    # Temporary code for testing:
    #    m = ModiaSimulationModel(model_name_of(instance), FDAE, x0, der_x0, jac;
    #                             maxSparsity=maxSparsity, nz=initial_m.nz_preInitial,
    #                             xNames=xNames, nc=sortedEquations.nc)
    #    m = ModiaSimulationModel(string(model_name_of(instance)), F, x0;
    #                             maxSparsity=maxSparsity, nc=1, nz=initial_m.nz_preInitial, jac=jac, x_fixed=diffstates)

    m = ModiaSimulationModel("test", ModuleFDAE.FDAE, x0) #, der_x0, nothing)
    #                        nc=nc)
    #=
    println("Call FDAE once")
    w = []
    Base.invokelatest(ModuleFDAE.FDAE, m, 0, x0, der_x0, r, w)
    println("FDAE called")
    @show r
    @show x0
    =#

    if length(x0) > 0
        println("\nSimulate")
#        t = linspace(0.0, 50, 1000)
        t = range(0.0, stop=50.0, length=1000)
        result = ModiaMath.ModiaToModiaMath.simulate(m, t; log=false, tolRel=1E-5)

#= Temporarily removed due to problem with PyPlot
        # Plot results
        t = result["time"]
        # figure()
        clf()
        # title(string(m))
        vars = keys(result)
        leg = String[]
        for v in vars
            if v != "time" && !(findfirst(isequal(v), "der") != notFound)
                plot(result["time"], result[v])
                push!(leg, v)
            end
        end
        legend(leg,  loc="center right")
        grid(true)
        xlabel("time [s]")
        # @show result vars leg
        # plot(result, Tuple(leg))
        #    (t_res,x_res,der_x_res) = ModiaMath.ModiaToModiaMath.getRawResult(m)
=#
        reuse=false
        for v in keys(result)
           if v != "time" && !(findfirst(isequal(v), "der") != notFound)
               ModiaMath.plot(result, v, reuse=reuse)
               reuse=true
           end
        end
    end
    # results = extract_results_ida(x_res, der_x_res, states, state_offsets, params)
    nothing
end


end