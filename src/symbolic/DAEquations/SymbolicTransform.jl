"""
Module for symbolic transformation of models.

* Developer: Hilding Elmqvist, Mogram AB  
* First version: July-September 2016
* License: MIT (expat)

Reference:
For matrix formulas, see:
Petersen, K.B, Pedersen M.S.: The Matrix Cookbook
https://www.math.uwaterloo.ca/~hwolkowi/matrixcookbook.pdf

"""
module SymbolicTransform

#using Debug
using ..Modia
using ..Instantiation
import ..Instantiation: GetField, This, Der, Symbolic, time_global, simulationModel_symbol
using Base.Meta: quot, isexpr
using DataStructures
using Unitful
@static if !(VERSION < v"0.7.0-DEV.2005")
    import LinearAlgebra
end
@static if VERSION < v"0.7.0-DEV.2005"
    notFound = 0
else
    notFound = nothing
end


using ..Synchronous
using ..BLTandPantelidesUtilities
#using ..Utilities
import ModiaMath
using ..ModiaLogging

export initSynchronousCounters, substituteClocked, solve, solveEquations, solveSortedEquations, differentiate, differentiateEquations, differentiateSortedEquations, residue_der, newDifferentiateEquations
export setDerAsFunction, setTimeInvariants, setRealStates

const logSolve = false
const logDifferentiate = true
const logDifferentiateVariable = false
const logSolved = false
const noUnits = true
const pushResiduals = false
const checkUnits = false

#=
@show logSolve
@show logSolved
=#

global derAsFunction = true
setDerAsFunction(v) = (global derAsFunction = v)

global nClocks
global nsamples
global nprevious
global npositive

function initSynchronousCounters()
    global nClocks = 0
    global nsamples = 0
    global nprevious = 0
    global npositive = 0
end

substituteClocked(previousVariables, ex) = ex

function substituteClocked(previousVariables, ex::Expr)
    global npositive
    if isexpr(ex, :quote)
        ex
    elseif isexpr(ex, :call) && ex.args[1] == Synchronous.Clock
        global nClocks += 1
        Expr(ex.head, ex.args[1], ex.args[2], simulationModel_symbol, nClocks)
    elseif isexpr(ex, :call) && ex.args[1] == Synchronous.sample
        global nsamples += 1
        Expr(ex.head, ex.args[1], ex.args[2], substituteClocked(previousVariables, ex.args[3]), simulationModel_symbol, nsamples)
    elseif isexpr(ex, :call) && ex.args[1] == Synchronous.previous
        global nprevious += 1
        push!(previousVariables, ex.args[2])
        Expr(ex.head, ex.args[1], ex.args[2], substituteClocked(previousVariables, ex.args[3]), simulationModel_symbol, nprevious)
    elseif isexpr(ex, :call) && ex.args[1] == Synchronous.positive
        global npositive += 1
        Expr(ex.head, ModiaMath.ModiaToModiaMath.positive!, ex.args[2], simulationModel_symbol, npositive)
    else    
        Expr(ex.head, [substituteClocked(previousVariables, arg) for arg in ex.args]...)
    end
end

zero = Float64(0.0)
one = Float64(1.0)

function add(e1, e2)
    if e1 == zero
        return e2
    elseif e2 == zero 
        return e1
    else
        Expr(:call, +, e1, e2)
    end
end

function sub(e1, e2)
    if typeof(e1) in [Float64, Int64] && typeof(e2) in [Float64, Int64]
        e1 - e2
    elseif e1 == zero
        return Expr(:call, -, e2)
    elseif e2 == zero 
        return e1
    else
        Expr(:call, -, e1, e2)
    end
end

function mult(e1, e2)
    if e1 == zero || e2 == zero  
        return zero
    elseif e1 == one
        return e2
    elseif e2 == one 
        return e1  
    else
        Expr(:call, *, e1, e2)
    end
end

function power(e1, e2)
    if e2 == one 
        return e1
    else
        Expr(:call, ^, e1, e2)
    end
end

# operators = [call, +, :+, -, :-, *, :*, (==), ^, :^, :(=), \, /, :/, transpose]
operators = [:call, +, :+, -, :-, *, :*, (==), ^, :^, :(=), \, /, :/, transpose] # v 0.6

checkincidence(s::Function) = s in operators ? nothing : s 
checkincidence(e::Instantiation.GetField) = e.name
# checkincidence(e::Instantiation.Der) = e.base
# checkincidence(e::Expr) = map(checkincidence, e.args)
checkincidence(e::Number) = nothing
checkincidence(e::LineNumberNode) = nothing
checkincidence(x) = x

# From: https://rosettacode.org/wiki/Flatten_a_list#Julia
@static if VERSION < v"0.7.0-DEV.2005"
    flat(A) = mapreduce(x -> isa(x, Array) ? flat(x) : x, vcat, [], A)
else
    flat(A) = mapreduce(x -> isa(x, Array) ? flat(x) : x, vcat, A, init=[])
end

function getincidence(e)
    # @show e typeof(e)
    if isa(e, Function)
        incidence = [e in operators ? nothing : e] 
    elseif typeof(e) <: Unitful.Quantity
        # dump(e)
        incidence = nothing
    elseif typeof(e) == Instantiation.GetField
        incidence = [e, e.name]
    elseif typeof(e) == Instantiation.Der
        incidence = e # getincidence(e.base)
    elseif typeof(e) == Expr
        if e.head == :call && e.args[1] == :der  
            incidence = e  
        elseif !derAsFunction && e.head == :call && e.args[1] == Synchronous.previous
            incidence = map(getincidence, e.args[3:end])
        else
            incidence = map(getincidence, e.args[2:end])
        end
    elseif e == time_global || e == :time
        incidence = nothing
    elseif typeof(e) == Symbol
        # if e in [:kg, :m, :s, :N]
            # incidence = nothing
        # else
        incidence = e 
        # end
    elseif typeof(e) in [Number, Float64, Int64, LineNumberNode]
        incidence = nothing
    # elseif typeof(e) in [block]
            # incidence = checkincidence(e.args[2])
    elseif e != nothing
        incidence = map(checkincidence, e)
    else
        incidence = nothing
    end

    if isa(incidence, Array)
      ### This allocates a lot of memory and should be rewritten.
        incidence = flat(incidence)
    else
        incidence = [incidence]
    end
    incidence = filter(e -> e != nothing, incidence)
    # println("getincidence: ", incidence)
    return incidence
end


find_incidence!(incidence::Array{Any,1}, ex) = nothing
find_incidence!(incidence::Array{Any,1}, der::Der) = (push!(incidence, der); nothing)
find_incidence!(incidence::Array{Any,1}, get::GetField) = (push!(incidence, get.name); nothing)
  
function find_incidence!(incidence::Array{Any,1}, ex::Expr)
    if !isexpr(ex, :quote)
        if ex.head == :call
            if ex.args[1] in [(Synchronous.previous), :(Synchronous.previous)]
                for arg in ex.args[3:end]    # 3:end ?
                    find_incidence!(incidence, arg)
                end   
            elseif ex.args[1] == :der  
                println(":der incidence")
                push!(incidence, ex)  
            end
        else
            for arg in ex.args
                find_incidence!(incidence, arg)
            end
        end
    end
    nothing
end

#=
inc = []
find_incidence!(inc, LHS)
if testIncidence(inc, x)
=#

function testIncidence(e, x)
    inc = getincidence(e)
    return findfirst(isequal(x), inc) != notFound
end
  

# --------------------------------------------------------------

function solve(eq::Expr, x)
    #=
    Solve equation eq for x 

    Rewriting rules:
    ex: expression containing x
    e = ex  ->  ex = e
    -ex = e  ->  ex = -e

    ex + e1 = e2  ->  ex = e2 - e1
    e1 + ex = e2  ->  ex = e2 - e1

    ex - e1 = e2  ->  ex = e2 + e1
    e1 - ex = e2  ->  ex = e1 - e2

    ex * e1 = e2  ->  ex = e2 / e1
    e1 * ex = e2  ->  ex = e1 \ e2

    ex / e1 = e2  ->  ex = e2 * e1
    e1 / ex = e2  ->  e1 = e2 * ex  ->  ex = e2 \ e1

    ex \ e1 = e2  ->  e1 = ex * e2  ->  ex = e1 / e2
    e1 \ ex = e2  ->  ex = e1 * e2
 
    transpose(ex) = e  ->  ex = transpose(e)
    =#

    #  @assert eq.head in [:(=), :(:=)] "solve only solves equations"
    if !(eq.head in [:(=), :(:=)])
        return (eq, false)
    end
    
    if logSolve
        println("\nSOLVE: ", string(x), " from:   ", prettyPrint(eq)) 
    end
    
    LHS = eq.args[1]
    RHS = eq.args[2]

    if testIncidence(LHS, x) && testIncidence(RHS, x)
        println("Multiple occurances of unknown: ", prettyPrint(eq))
        return eq, false
    elseif RHS == x || testIncidence(RHS, x)
        # Unknown x in RHS
        if logSolve println("  right incidence, swap") end
        # Swap LHS and RHS
        eq, solved = solve(Expr(:(=), RHS, LHS), x)
        return eq, solved 
    elseif LHS == x || testIncidence(LHS, x)
        # Unknown x in LHS
        if LHS == x || typeof(LHS) == Symbol || typeof(LHS) == Instantiation.GetField || typeof(LHS) == Instantiation.Der
            if logSolve println("  Already solved: ", prettyPrint(eq)) end
            return eq, true
        elseif LHS.head == :block && (isa(LHS.args[1], LineNumberNode) || isexpr(LHS.args[1], :line))
            eq = Expr(:(=), LHS.args[2], RHS)
            eq, solved = solve(eq, x)
            return eq, solved
        end

        op = LHS.args[1]
        arguments = LHS.args[2:end]
        er = RHS 

        if size(arguments, 1) == 1 && op in [-, :-, :transpose]
            # Unary operators - and transpose
            #  - ex = er  ->  ex = -er
            #  transpose(ex) = er  ->  ex = transpose(er)
            ex = LHS.args[2]
            eq = Expr(:(=), ex, Expr(:call, op, er))       

            if ex != x
                eq, solved = solve(eq, x)
                return eq, solved
            else
                return eq, true
                dump(ex)
                dump(x)
                error("Should not happen")
            end

        elseif op in [*, :*, +, :+, -, :-, /, :/]   # ^, :^]
            e1 = arguments[1]
            if size(arguments, 1) == 2
                e2 = arguments[2]
            else
                # e2 = Expr(:call, op, arguments[2:end])  # Not correct
                e2 = Expr(:call) # , op, arguments[2])
                e2.args = [op; arguments[2:end]]  # Fix
            end
      
            if testIncidence(e1, x) && testIncidence(e2, x)
                println("Multiple occurances of unknown: ", prettyPrint(LHS))
                return eq, false
            elseif testIncidence(e1, x)
                if logSolve println("  left expression incidence") end
                ex = e1
                if op in [*, :*] 
                    # ex * e2 * ... = er  ->  ex = er / (e2 * ...)
                    eq = Expr(:(=), ex, Expr(:call, :/, er, e2))
                elseif op in [+, :+] 
                    # ex + e2 + ... = er  ->  ex = er - (e2 + ...)
                    eq = Expr(:(=), ex, sub(er, e2))
                elseif op in [-, :-] 
                    # ex - e2 - ... = er  ->  ex = er + e2 + ...
                    if size(arguments, 1) > 2
                        e2.args[1] = +
                    end
                    eq = Expr(:(=), ex, add(er, e2))
                elseif op in [/, :/] 
                    # ex / e1 / ... = er  ->  ex = er * e1 * ...
                    if size(arguments, 1) > 2
                        e2.args[1] = *
                    end
                    eq = Expr(:(=), ex, mult(er, e2))
                elseif op in [^, :^] 
                    # ex ^ 2 = er  ->  ex = sqrt(er)  # Generalize, could be negative
                    # eq = Expr(:(=), ex, Expr(:call, :-, Expr(:call, :sqrt, er)))
                    return eq, false
                end
            else
                if logSolve println("  right expression incidence") end
                ex = e2
                if op in [*, :*]
                    # e1 * ex * ... = er  ->  ex * ... = e1 \ er
                    eq = Expr(:(=), ex, Expr(:call, :\, e1, er))
                elseif op in [+, :+]
                    # e1 + ex + ... = er  ->  ex + ... = er - e1
                    eq = Expr(:(=), ex, sub(er, e1))
                elseif op in [-, :-]
                    # e1 - ex - ... = er  ->  ex + ... = e1 - er
                    if size(arguments, 1) > 2
                        ex.args[1] = +
                    end
                    eq = Expr(:(=), ex, sub(e1, er))
                elseif op in [/, :/]
                    # e1 / ex / ... = er  ->  ex * ... = e1 / er
                    if size(arguments, 1) > 2
                        ex.args[1] = *
                    end
                    eq = Expr(:(=), ex, Expr(:call, /,  e1, er))
                end
            end
    
            if ex != x
                eq, solved = solve(eq, x)
                return eq, solved
            else
                if logSolve println("  Solved: ", prettyPrint(eq)) end
                return eq, true
            end
        else
            if logSolve 
                println("Equations involving this operator are not solved: ", op)
                println("  Not solved: ", prettyPrint(eq))
            end
            # dump(eq)
            return eq, false
        end

        if logSolve println("  Never solved: ", prettyPrint(eq)) end
        error("Should not happen")
        return eq, false
    else
        # println("Should not happen:")
        # dump(eq.args[1])
        # dump(eq.args[2])
        # dump(x)
        return eq, false
        error("Should not happen")  # Check up!!!
    end
  
    if logSolve println("  Not solved: ", prettyPrint(eq)) end
    error("Should not happen")
    return eq, false
end

"""
    solveEquations(equations, variables, indices, assign, A, B)

Solve assigned equations.

* `equations`: infix string for original equations
* `variables`: infix string for original variables
* `indices`: indices for the equations to be printed
* `assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `A`: A[j] = if V[k] = der(V[j]) then k else 0
* `B`: B[i] = if E[l] = der(E[l]) then l else 0

"""
function solveEquations(equations, variables, indices, assign, A, B)
    (orgIndexVar, derOrderVar) = invertDer(A)
    (orgIndexEqu, derOrderEqu) = invertDer(B)
    (assignedVar,) = invertAssign(assign)
    
    for i in indices
        j = assignedVar[i]
        if j > 0    
            if derOrderVar[j] == 1
                print("der(")
            elseif derOrderVar[j] > 1
                print("der", derOrderVar[j], "(")
            end
    
            print(variables[orgIndexVar[j]])
    
            if derOrderVar[j] > 0
                print(")")
            end
        end
        print(": ")
      
        if derOrderEqu[i] == 1
            print("DER() ")
        elseif derOrderEqu[i] > 1
            print("DER", derOrderEqu[i], "() ")
        end
    
        println(equations[orgIndexEqu[i]])
        println(solve(equations[i], variables[j]))
    end
end


"""
    solveSortedEquations(equations, variables, components, assign, A, B)

Solve sorted equations.

* `equations`: infix string for original equations
* `variables`: infix string for original variables
* `components`: cell array of components. Each component is a list of indices to E-nodes
* `assign`: assign[j] contains the E-node to which V-node j is assigned or 0 if V-node j not assigned
* `A`: A[j] = if V[k] = der(V[j]) then k else 0
* `B`: B[i] = if E[l] = der(E[l]) then l else 0

"""
function solveSortedEquations(equations, variables, components, assign, A, B)
    println("[assigned variable]: [differentiation] equation")
    println("Strongly connected components are enclosed in []")
    
    for c in components
        if length(c) > 1
            println("[")
        end
    
        solveEquations(equations, variables, c, assign, A, B)
    
        if length(c) > 1
            println("]")
        end
    end
end

# ----------------------------------------------------------------------------------------

# Derivatives of functions
divide_der_1(x, y) = 1 / y
divide_der_2(x, y) = x / y^2

power_der_1(x, y) = y * x^(y - 1)
power_der_2(x, y) =  x^y * log(x)
  
exp_der(x) = exp(x)
log_der(x) = 1 / x
sin_der(x) = cos(x)
cos_der(x) = -sin(x)
tan_der(x) = 1 / cos(x)^2
arcsin_der(x) = 1 / sqrt(1 - x^2)
arccos_der(x) = -1 / sqrt(1 - x^2)
arctan_der(x) = 1 / (1 + x^2)
sqrt_der(x) = 0.5 / sqrt(x)
sqrt_der_der(x) = -0.25 / (x * sqrt(x))
 

global timeInvariants = []
setTimeInvariants(tinv) = (global timeInvariants; timeInvariants = tinv)

global realStates = []
setRealStates(st) = (global realStates; realStates = st)
  
global dummyDerivatives
  
differentiate(e::Irrational) = zero
  
function differentiate(e)
    if logDifferentiate
        # print("\nDIFFERENTIATE: "); println(prettyPrint(e))
        # @show typeof(e)
    end
    
    if typeof(e) in [Float64, Int64, String, Bool]
        # diff = if ! derAsFunction; :($zero*$e) else zero*e end
        diff = if !derAsFunction; :($zero) else zero * e end
        if !noUnits 
            diff = diff / SIUnits.Second
        end
    elseif typeof(e) <: Unitful.Quantity ||  typeof(e) <: Unitful.Unitlike
        diff = zero * e / Unitful.s
    elseif typeof(e) in [Array{Float64}, Array{Int64}]
        # diff = :($zero*$e) # zero*e
        diff = :($zero) # zero*e
    elseif e == :time
        diff = one
    elseif e in timeInvariants
        # println("TIME INVARIANT" )
        diff = zero
    elseif typeof(e) == Symbol
        diff = Expr(:call, :der, e)
    elseif typeof(e) == Instantiation.GetField
        # @show e.name
        if e.name in timeInvariants
            # println("TIME INVARIANT" )
            # diff = :($zero*$e)
            diff = :($zero)
        else
            if logDifferentiateVariable
                      print("DIFFERENTIATE VARIABLE:    "); println(e)
            end
            # println("\nGetField:")
            # @show e
            # dump(e)
            i = findfirst(isequal(e), realStates)
            
            if i != notFound
                diff = Der(e) # / Unitful.s 
                # @show realStates i e diff
            else
                # diff = GetField(This(), Symbol("DER("*string(e.name)*")"))
                # dummyDerivatives[Symbol("DER("*string(e.name)*")")] = nothing
                diff = GetField(This(), Symbol("der_" * string(e.name))) # / Unitful.s # With units
                dummyDerivatives[Symbol("der_" * string(e.name))] = nothing
                
                if logDifferentiateVariable
                    println("  Dummy derivative: ", diff)
                end
            end
            
            # diff = Der(e)
            if logDifferentiateVariable
                print("  Differentiated variable: "); println(diff)
            end
        end

    elseif typeof(e) == Instantiation.Der
        if logDifferentiateVariable
                print("DIFFERENTIATE DERIVATIVE:    "); println(e)
        end

#        @show e e.base e.base.name realStates
        state = GetField(This(), Symbol("der_" * string(e.base.name)))
#        @show state        
        i = findfirst(isequal(state), realStates)
#        @show i
        if i == notFound
            diff = GetField(This(), Symbol("der_der_" * string(e.base.name))) 
#            @show diff
            if !noUnits 
                diff = diff / Unitful.s
            end
            dummyDerivatives[Symbol("der_der_" * string(e.base.name))] = nothing
        else
            diff = Der(state)
            dummyDerivatives[Symbol("der_" * string(e.base.name))] = nothing
        end
        
        # @show e diff
        if logDifferentiateVariable
            print("  Differentiated derivative: "); println(diff)
        end

    elseif e.head in [:(=), :(:=)] # Equation
        diff = Expr(e.head, differentiate(e.args[1]), differentiate(e.args[2])) 
        if logDifferentiate
            logModia("DIFFERENTIATE:    "); loglnModia(prettyPrint(e))
            logModia("  Differentiated: "); loglnModia(prettyPrint(diff)); loglnModia()
        end
    elseif e.head == :block && (isa(e.args[1], LineNumberNode) || isexpr(e.args[1], :line))
        diff = differentiate(e.args[2])
    elseif e.head in [:vect, Symbol("vect"), :row, Symbol("row"), :vcat, Symbol("vcat")]
        diff = Expr(e.head)
        for i in 1:length(e.args)
            push!(diff.args, differentiate(e.args[i]))
        end
        # @show e
        # @show diff
    elseif e.head in [:ref]
        diff = copy(e)
        diff.args[1] = differentiate(e.args[1])  # Probably not complete.
    elseif e.head in [:comprehension]
        diff = copy(e)
        diff.args[2].args[2] = differentiate(e.args[2].args[2])
    elseif e.head in [:if]
        diff = Expr(e.head, e.args[1], differentiate(e.args[2]), differentiate(e.args[3]))          
    elseif e.head == :call
        op = e.args[1]
        arguments = e.args[2:end]
        if length(arguments) == 1 && op in [:der, :transpose, transpose, Symbol("'"), :residue, Modia.BasicStructuralTransform.residue, :skew, Modia.BasicStructuralTransform.skew] # skewCoords
            diffarg = differentiate(arguments[1])
            # @show e
            # @show diffarg
            if diffarg == zero
                diff = zero
            else
                diff = Expr(:call, op, diffarg)  
            end
            # @show diff
        elseif (VERSION < v"0.7.0-DEV.2005" && op in [:zeros, zeros, :eye, eye] ) ||
              ( ! (VERSION < v"0.7.0-DEV.2005") && op in [:zeros, zeros, :eye, LinearAlgebra.I, :identity, identity] )
            diff = zero # !!!!
            # diff = Expr(:call, *, e, zero)  # = e*zero
        elseif op in [+, :+, -, :-]
            # der(e1 + e2 + e3) = der(e1) + der(e2) + der(e3)
            # der(e1 - e2 - e3) = der(e1) - der(e2) - der(e3)
            # diff = Expr(:call, op, [differentiate(arg) for arg in arguments])  # Not correct
            # diff = Expr(:call, op)
            # diff.args = [op; [differentiate(arg) for arg in arguments]]  # Fix   
            diff = Expr(:call, op)
            firstNonZero = false
            
            for i in 1:length(arguments)
                term = differentiate(arguments[i])
                if term != zero 
                    if i == 1
                        firstNonZero = true
                    end
                    push!(diff.args, term)
                end
            end
            
            if length(diff.args) == 1
                diff = zero
            elseif op in [+, :+] && length(diff.args) == 2
                # The cases:
                # der(x+5) -> der(x)        
                diff = diff.args[2]
            elseif op in [-, :-] && length(arguments) > 1 && length(diff.args) == 2 && firstNonZero
                # The cases:
                # der(x-5) -> der(x)
                # der(5-x) -> -der(x)
                # der(-x) -> -der(x)
                diff = diff.args[2]
            end
        elseif op in [*, :*]
            # der(e1 * e2 * e3) = der(e1)*e2*e3 + e1*der(e2)*e3 + e1*e2*der(e3)
            arguments = e.args[2:end]
            diff = Expr(:call, +)
            for i in 1:length(arguments)
                terms = []
                for j in 1:length(arguments)
                    term = if i == j; differentiate(arguments[j]) else arguments[j] end
                    if term == zero  
                        terms = []
                        break
                    else
                        push!(terms, term)
                    end
                end
            
                if terms != []
                    product = Expr(:call, *)
                    for t in terms
                        push!(product.args, t)
                    end
            
                    if length(diff.args) > 0
                        push!(diff.args, product)   
                    end             
                end
            end   
            
            if length(diff.args) == 1 # Only + operator
                diff = zero
            elseif length(diff.args) == 2 # Only + zero, remove +
                diff = diff.args[2]
            end
        elseif op in [/, :/]
            # der(e1/e2) = der(e1)/e2 + e1/e2^2*der(e2)
            e1 = arguments[1]
            e2 = arguments[2]
            d1 = differentiate(e1)
            d2 = differentiate(e2)
            
            if d1 == zero
                diff = :($e1 / $e2^2 * $d2)
            elseif d2 == zero
                diff = :($d1 / $e2)
            else
                diff = :($d1 / $e2 + $e1 / $e2^2 * $d2)
            end
        elseif op in [^, :^] && typeof(arguments[2]) in [Int64]
            # der(e1 ^ n) = n*e1^(n-1)*der(e1)
            e1 = arguments[1]
            n = arguments[2]
            s = sub(n, 1)
            p = power(e1, s)
            d = differentiate(e1)
            diff = :($n * $p * $d)
        elseif op in [^, :^]
            # der(e1^e2) = e2*e1^(e2-1)*der(e1) + e1^e2*log(e1)*der(e2)
            e1 = arguments[1]
            e2 = arguments[2]
            d1 = differentiate(e1)
            d2 = differentiate(e2)
            if d1 == zero
                diff = :($e1^$e2 * log($e1) * $d2)
            elseif d2 == zero
                diff = :($e2 * $e1^($e2 - 1) * $d1)
            else        
                diff = :($e2 * $e1^($e2 - 1) * $d1 + $e1^$e2 * log($e1) * $d2)
            end
        elseif op in [log, :log]
            # der(log(x)) = 1/x*der(x)
            e1 = arguments[1]
            d = differentiate(e1)
            diff = :(1 / $e1 * $d)
        elseif op in [exp, :exp]
            # der(exp(x)) = exp(x)*der(x)
            e1 = arguments[1]
            d = differentiate(e1)
            diff = :(exp($e1) * $d)
        elseif op in [cos, :cos]
            # der(cos(e1)) = -sin(e1)*der(e1)
            e1 = arguments[1]
            d = differentiate(e1)
            diff = :(-sin($e1) * $d)
        elseif op in [sin, :sin]
            # der(sin(e1)) = cos(e1)*der(e1)
            e1 = arguments[1]
            d = differentiate(e1)
            diff = :(cos($e1) * $d)
        elseif op in [tan, :tan]
            # der(tan(x)) = 1/cos(x)^2*der(x)
            e1 = arguments[1]
            d = differentiate(e1)
            diff = :(1 / cos($e1)^2 * $d)
        elseif op in [asin, :asin]
            # der(arcsin(x)) = 1/sqrt(1-x^2)*der(x)
            e1 = arguments[1]
            d = differentiate(e1)
            diff = :(1 / sqrt(1 - $e1^2) * $d)
        elseif op in [acos, :acos]
            # der(arccos(x)) = -1/sqrt(1-e1^2)*der(x)
            e1 = arguments[1]
            d = differentiate(e1)
            diff = :(-1 / sqrt(1 - $e1^2) * $d)        
        elseif op in [atan, :atan]
            # der(arctan(x)) = 1/(1+x^2)*der(x)
            e1 = arguments[1]
            d = differentiate(e1)
            diff = :(1 / (1 + $e1^2) * $d)
        else
            # der(f(e1, e2, e3)) = f_der_1(e1, e2, e3)*der(e1) + f_der_2(e1, e2, e3)*der(e2) + f_der_3(e1, e2, e3)*der(e3)
            arguments = e.args[2:end]
            diff = Expr(:call, +)
            if op in [/, :/]
                op = :divide
            elseif op in [^, :^]
                op = :power
            end
            
            for i in 1:length(arguments)
                arg_der = differentiate(arguments[i])
                if arg_der != zero
                    path = split(string(op), ".")
                    if length(path) >= 2
                        mod = path[end-1]
                    else 
                        mod = ""
                    end
                    func = path[end]
                    f_der = Symbol(string(func * "_der_", i))
                    
                    if mod == ""
                        println("Derivative function ", string(f_der), " not found.")
                    elseif ! (f_der in names(getfield(Main, Symbol(mod))))
                        error("Derivative function ", string(f_der), " not found.")
                    else
                        f_der = getfield(getfield(Main, Symbol(mod)), f_der)
                    end
                    f_der_i = Expr(:call, f_der)
                    for a in arguments
                        push!(f_der_i.args, a)
                    end
                
                    term = mult(f_der_i, arg_der) 
                    push!(diff.args, term)      
                end
            end       

            if length(diff.args) == 1
                push!(diff.args, zero)
            elseif length(diff.args) == 2
                diff = diff.args[2]
            end
        end
    else
        println("Error: ", prettyPrint(e))
        println("Not able to differentiate the above expression yet. ")
        dump(e)
        diff = nothing
    end
    # println("diff=", diff)
    return diff
end


function differentiateEquations(equations, variables, indices, assign, A, B, residualIndex)
    (orgIndexVar, derOrderVar) = invertDer(A)
    (orgIndexEqu, derOrderEqu) = invertDer(B)
    (assignedVar,) = invertAssign(assign)
    nSolved = 0
    component = []
    algebraic = []
    for i in indices
        if logSolve
            println("\nNEXT EQUATION")
        end
     
        j = assignedVar[i]
        assigned = j > 0
        if assigned
            if logSolved
                if derOrderVar[j] == 1
                    print("der(")
                elseif derOrderVar[j] > 1
                    print("der", derOrderVar[j], "(")
                end
     
                print(variables[orgIndexVar[j]])
                if derOrderVar[j] > 0
                    print(")")
                end
            end 
        end

        if logSolved
            print(": ")
        
            if derOrderEqu[i] == 1
                print("DER() ")
            elseif derOrderEqu[i] > 1
                print("DER", derOrderEqu[i], "() ")
            end
            println(prettyPrint(equations[orgIndexEqu[i]]))
        end
    
        diff = equations[orgIndexEqu[i]]
        for d in 1:derOrderEqu[i]
            diff = differentiate(diff)
        end
     
        if logSolved
            println(prettyPrint(diff))
        end
    
        if assigned
            solveFor = variables[orgIndexVar[j]]
            solveForDummyDer = solveFor
            stateVar = GetField(This(), solveFor)
            # @show stateVar
            for d in 1:derOrderVar[j]
                # @show d
                if derAsFunction
                    solveFor = Expr(:call, :der, solveFor)
                    solveForDummyDer = Symbol("der_" * string(solveForDummyDer))
                else
                    # @show solveFor
                    # solveFor = Der(GetField(This(), solveFor))
                    # @show solveFor
                    # k = findfirst(realStates, stateVar)
                    # @show k
                    # if k == 0
                    if !(stateVar in realStates) ### Patch
                        # @show solveForDummyDer
                        # @show realStates
                        if d == 1
                            # solveForDummyDer = GetField(This(), Symbol("der_"*string(solveForDummyDer)))
                            solveForDummyDer = GetField(This(), Symbol("der_" * string(solveFor)))
                        else
                            solveForDummyDer = GetField(This(), Symbol("der_der_" * string(solveFor)))            
                        end

                        if logDifferentiateVariable
                            println("DUMMY DERIVATIVE: ", solveForDummyDer)
                        end
                    else
                        solveForDummyDer = if !Modia.BasicStructuralTransform.useKinsol; nothing else Der(GetField(This(), solveFor)) end ### nothing # solveFor
                    end
                end
            end
            solveFor = solveForDummyDer
        else
            solveFor = nothing
        end

        if assigned && solveFor != nothing && length(indices) == 1
            (solution, solved) = solve(diff, solveFor)
        else
            solution = diff
            # solved = diff.head == :(:=)  # false  # To allow torn equations
            solved = false
        end
      
        if (assigned && solved && length(indices) == 1)  ||  !derAsFunction
      
    else
            push!(algebraic, solveFor)
            residualIndex += 1
            if noUnits && pushResiduals
                push!(component, (:(show = @show length(RESIDUALS)), true))
                push!(component, (:(show = @show $(diff.args[1]) - $(diff.args[2])), true))
                solution = :(dummy = push!(RESIDUALS, $(diff.args[1]) - $(diff.args[2])...))
            elseif noUnits
                solution = :(RESIDUALS[$residualIndex] = $(diff.args[1]) - $(diff.args[2]))
            else
                solution = :(RESIDUALS[$residualIndex] = float($(diff.args[1]) - $(diff.args[2])))
            end
            solved = false
        end
    
        push!(component, (solution, solved))
        if solved
            nSolved += 1
        end

        if logSolved    
            if !solved 
                print("NOT SOLVED: ")
                println()
                @show solveFor assigned solution
            end
            println(prettyPrint(solution))
        end
    end
    return component, algebraic, residualIndex
end

function differentiateSortedEquations(equations, variables, components, assign, A, B)
    if logSolved    
        println("[assigned variable]: [differentiation] equation")
        println("Strongly connected components are enclosed in []")
    end

    global dummyDerivatives = OrderedDict()
    solvedComponents = []
    algebraic = []
    residualIndex = 0
    
    for c in components
        if logSolved && length(c) > 1
            println("[")
        end
    
        (solvedComponent, newAlgebraic, residualIndex) = differentiateEquations(equations, variables, c, assign, A, B, residualIndex)
        push!(solvedComponents, solvedComponent)
    
        for a in newAlgebraic
            push!(algebraic, a)
        end
  
        if logSolved && length(c) > 1
            println("]")
        end
    end
    
    # @show algebraic
    # @show solvedComponents
    # println("End of symbolic processing.")
    return solvedComponents, algebraic, residualIndex
end

function calcLinearCoeffcients(e)
  
end

function getLinearCoeffcients(e)
#=
  
=#
    linearCoefficients = Dict()
    vars = []
    BasicStructuralTransform.findIncidence!(vars, :(5x + 5 * (2y + x * 8 + z / 3 = 0)))
    linearCoefficients[:1] = 0
    
    for i in vars
        linearCoefficients[i] = 0
    end
    @show linearCoefficients
end

function newDifferentiateEquations(equations, variables, A, B, ESorted, ESolved)
    global dummyDerivatives = OrderedDict()
    println("newDifferentiateEquations")
    (orgIndexVar, derOrderVar) = invertDer(A)
    (orgIndexEqu, derOrderEqu) = invertDer(B)
    nSolved = 0
    component = []
    algebraic = []
    
    for k in 1:length(ESorted)
        i = ESorted[k]
        if logSolve
            println("\nNEXT EQUATION: $i")
        end
        j = ESolved[k]
        # @show ESolved[k]
        assigned = j > 0
        if assigned
            if logSolved
                if derOrderVar[j] == 1
                    print("der(")
                elseif derOrderVar[j] > 1
                    print("der", derOrderVar[j], "(")
                end
        
                print(variables[orgIndexVar[j]])
        
                if derOrderVar[j] > 0
                    print(")")
                end
            end 
        end

        if logSolved
            print(": ")
        
            if derOrderEqu[i] == 1
                print("DER() ")
            elseif derOrderEqu[i] > 1
                print("DER", derOrderEqu[i], "() ")
            end
            println(prettyPrint(equations[orgIndexEqu[i]]))
        end
    
        diff = equations[orgIndexEqu[i]]

        for d in 1:derOrderEqu[i]
            diff = differentiate(diff)
        end
        
        if logSolved
            println(prettyPrint(diff))
        end
        
        # @show assigned
        if assigned
            solveForOrg = variables[orgIndexVar[j]]
            solveFor = variables[orgIndexVar[j]]

            for d in 1:derOrderVar[j]
                # @show d
                # @show solveFor
                # dump(solveFor)

                if d == 1
                    # solveForDummyDer = GetField(This(), Symbol("der_"*string(solveForDummyDer)))
                    if false ###
                        solveFor = GetField(This(), Symbol("der_" * string(solveFor)))
                    else
                        solveFor = Der(GetField(This(), solveFor))
                    end          
                elseif d == 2
                    solveFor = GetField(This(), Symbol("der_der_" * string(solveForOrg)))            
                else
                    solveFor = GetField(This(), Symbol("der_der_" * string(solveFor)))            
                end
            end
 
        else
            solveFor = nothing
        end
    
        # @show solveFor
        if assigned && solveFor != nothing 
            # @show solveFor
            (solution, solved) = solve(diff, solveFor)
            # @show solved
            if !solved
                # dump(solveFor)
                # dump(diff)
            end
        else
            solution = diff
            # solved = diff.head == :(:=)  # false  # To allow torn equations
            solved = false
        end
      
        if (assigned && solved) # ||  ! derAsFunction
      
        else
            push!(algebraic, solveFor)
            residualIndex = -j
        
            if noUnits && pushResiduals
                push!(component, (:(show = @show length(_r)), true))
                push!(component, (:(show = @show $(diff.args[1]) - $(diff.args[2])), true))
                solution = :(dummy = push!(_r, $(diff.args[1]) - $(diff.args[2])...))
            elseif noUnits
                solution = :(_r[$residualIndex] = $(diff.args[1]) - $(diff.args[2]))
            else
                solution = :(_r[$residualIndex] = float($(diff.args[1]) - $(diff.args[2])))
            end
            solved = false
        end
    
        push!(component, (solution, solved))
        
        if solved
            nSolved += 1
        end
        
        if logSolved    
            if !solved 
                print("NOT SOLVED: ")
            end
            println(prettyPrint(solution))
        end
    end
  
    return component
end

end
  
