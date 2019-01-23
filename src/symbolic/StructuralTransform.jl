"""
Module for structural transformation of models including alias handling and size and type deduction.

* Developer: Hilding Elmqvist, Mogram AB  
* First version: July-August 2016
* Copyright (c) 2016-2019: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)
"""
module StructuralTransform

using ..BLTandPantelides
using ..BLTandPantelidesUtilities
using ..Utilities
using ..BasicStructuralTransform
using ..Instantiation
import ..Instantiation: GetField, This, Der, Symbolic, time_global, get_dims, parameter
import ..Execution: get_value, Subs, subs, Variable, split_variables, vars_of, GetField
using Base.Meta: quot, isexpr
using DataStructures
using ..Synchronous
using ..SymbolicTransform
using ..ModiaLogging

@static if VERSION < v"0.7.0-DEV.2005"
    notFound = 0
else
    notFound = nothing
end

#export elaborate, prettyPrint, simulateModel, skewCoords, skew, residue, hide, modiaCross, showExpr, transformModel
export residue, residue_der, hide  # from BasicStructuralTransform to models

global aliasElimination
global deduceSizes
const log = true
const extendedLog = false
const debug = false

const tSize = Tuple{Int64,Vararg{Int64,N} where N}    
const tSizes = Array{Tuple{Int64,Vararg{Int64,N} where N},1}    

export setOptions
function setOptions(options) 
    global aliasElimination = false
    if haskey(options, :aliasElimination)
        global aliasElimination = options[:aliasElimination]
        @show aliasElimination
        delete!(options, :aliasElimination)
    end
    global deduceSizes = true
    if haskey(options, :deduceSizes)
        global deduceSizes = options[:deduceSizes]
        @show deduceSizes
        delete!(options, :deduceSizes)
    end
end

findStates!(states, deriv::Vector, ex) = nothing
findStates!(states, deriv::Vector, der::Der) = (push!(states, der.base); push!(deriv, der); nothing)
function findStates!(states, deriv::Vector, ex::Expr)
    if !isexpr(ex, :quote)
        for arg in ex.args
            findStates!(states, deriv, arg)
        end
    end
    nothing
end

findIncidence!(incidence::Array{Any,1}, ex) = nothing
findIncidence!(incidence::Array{Any,1}, der::Der) = (push!(incidence, der); nothing)
findIncidence!(incidence::Array{Any,1}, get::GetField) = (push!(incidence, get.name); nothing)
function findIncidence!(incidence::Array{Any,1}, ex::Expr)
    if !isexpr(ex, :quote)
        if ex.head == :call && ex.args[1] == hide
    elseif ex.head == :call && (ex.args[1] == Synchronous.previous || ex.args[1] == Synchronous.sample)
            for arg in ex.args[3:end]  # Skip first argument to previous() and sample()
                findIncidence!(incidence, arg)
            end        
        else
            for arg in ex.args
                findIncidence!(incidence, arg)
            end
        end
    end
    nothing
end
  
function findNonStateVariables(src::VariableDict)
    nonStateVariables = []
    for (name, var) in src
        if isa(var, Variable) && !var.state
            push!(nonStateVariables, name)
        end
    end
    return nonStateVariables
end


# Alias elimination

#=
The set of equations is searched for trivial equations of the form:
v = number  # replace v by number
v1 = v2     # replace v1 by v2
x = v       # replace v by x if der(x) appear but not der(v)

v1 := v2
v3 = v2 -> 
=#

function equivalence(pairs)
    equ = []
    for p in pairs
        for e in equ
            if intersect(p, e) != []
                if p[1] in e && p[2] in e
                    println("Redundant equation: $(p[1]) = $(p[2])")
                end
                union!(p, e)
                setdiff!(e, e)
            end
        end
        push!(equ, p)
    end
    equivalenceClasses = []
    for e in equ
        if e != []
            push!(equivalenceClasses, e)
        end
    end
    equivalenceClasses
end

function findAliases2!(nonAliasEquations, aliases, equations, variables, defined)
    equivalencePairs = []
    for eq in equations
        if !isexpr(eq, :quote)
            if eq.head in [:(=), :(:=)]
                e1 = eq.args[1]
                e2 = eq.args[2]
                if typeof(e1) == GetField && typeof(e2) == GetField 
                    push!(equivalencePairs, [e1.name, e2.name])
                else
                    push!(nonAliasEquations, eq)
                end
            end
        end
    end
    equivalenceClasses = equivalence(equivalencePairs)
    for equiv in equivalenceClasses
        nonalias = nothing
        for v in equiv
            if typeof(variables[v]) != Variable || variables[v].variability <= parameter
                nonalias = GetField(This(), v)
                break
            end
        end
        if nonalias == nothing
            nonalias = GetField(This(), equiv[1]) # No non-Variable or parameter to define as nonalias, pick first
        end
        for v in equiv
            if v != nonalias.name
                aliases[GetField(This(), v)] = nonalias
                if typeof(variables[v]) != Variable || variables[v].variability <= parameter
                    if "$(variables[v])" != "$(variables[nonalias.name])"
                        println("Non consistent constraint: $v = $(nonalias.name): $(variables[v]) != $(variables[nonalias.name])")
                    end
                end
            end
        end
    end
end
        

function findAliases!(nonAliasEquations, aliases, eq::Expr, unknowns, defined)
    if !isexpr(eq, :quote)
        if eq.head in [:(=), :(:=)]
            e1 = eq.args[1]
            e2 = eq.args[2]
            @show e1 e2
            if typeof(e1) == GetField && typeof(e2) == GetField && haskey(unknowns, e1.name) && haskey(unknowns, e2.name) # && unknowns[e1.name].state && unknowns[e2.name].state
                # nonAliases = collect(values(aliases))
                # if ! haskey(aliases, e1) && ! (e1 in nonAliases) && ! haskey(aliases, e2) && ! (e2 in nonAliases)
                # if ! haskey(aliases, e1) && ! (e1 in nonAliases) && ! haskey(aliases, e2) && ! (e2 in nonAliases)
                if !haskey(aliases, e1) && getNonAliasVariable(aliases, e2, unknowns)[1] != e1 && if haskey(unknowns, e2.name); unknowns[e2.name].state else false end 
                    aliases[e1] = e2
                    unknowns[e2.name].state = unknowns[e1.name].state
                    # @show e2.name, unknowns[e2.name].state
                    println()
                elseif !haskey(aliases, e2) && getNonAliasVariable(aliases, e1, unknowns)[1] != e2 && if haskey(unknowns, e1.name); unknowns[e1.name].state else false end
                    aliases[e2] = e1
                    unknowns[e1.name].state = unknowns[e2.name].state
                    @show e1.name, unknowns[e1.name].state
                else 
                    println("Circular aliases between: $e1 and $e2")
#                    error("Aborting")
                end
            elseif typeof(e1) == GetField && typeof(e2) == GetField && (haskey(unknowns, e1.name) || haskey(unknowns, e2.name))
                println("One is unknown:")
                println(eq)
                if haskey(unknowns, e1.name)
                    var = e1.name
                    val = e2.name
                else
                    var = e2.name
                    val = e1.name
                end
                @show var
                @show aliases
                if haskey(aliases, var)
                    var = aliases[var]
                    @show aliases[var]
                end
                @show defined
                if ! haskey(defined, var)
                    defined[var] = val
                    @show defined
                    println("Adding: ", eq)
                    push!(nonAliasEquations, eq)
                elseif defined[var] != val
                    println("Non consistent overdetermined variable:")
                    println("$var = $val")
                    println("$var = $defined[var]")
                else
                    println("Not added: ", eq)
                end                
            else
                println(eq)
                push!(nonAliasEquations, eq)
            end
        end
    end
    nothing
end

function getNonAliasVariable(aliases, ex, unknowns)
    # println("getNonAliasVariable")
    # @show aliases ex unknowns
    if haskey(aliases, ex)
        (e2, state) = getNonAliasVariable(aliases, aliases[ex], unknowns)
        # @show e2
        return e2, unknowns[ex.name].state && state
    else
        return ex, if haskey(unknowns, ex.name); unknowns[ex.name].state else false end
    end
end

const AliasSubs = Dict{Symbolic,Symbolic}
  
substiteAlias(aliases::AliasSubs, ex, unknowns) = get(aliases, ex, ex)

function substiteAlias(aliases::AliasSubs, ex::Symbolic, unknowns)
    if haskey(aliases, ex)
        getNonAliasVariable(aliases, aliases[ex], unknowns)[1]
    else
        ex
    end
end

function substiteAlias(aliases::AliasSubs, ex::Der, unknowns) 
    if haskey(aliases, ex.base)
        (s, state) = getNonAliasVariable(aliases, aliases[GetField(This(), ex.base.name)], unknowns)
        Der(GetField(This(), s.name))
    else
        ex
    end
end

function substiteAlias(aliases::AliasSubs, ex::Expr, unknowns)
    if isexpr(ex, :quote)
        ex
    else
        Expr(ex.head, [substiteAlias(aliases, arg, unknowns) for arg in ex.args]...)
    end
end

function performAliasElimination!(flat_model, unknowns, params, equations)
    loglnModia("\nALIAS ELIMINATION")
    aliases = AliasSubs()
    nonAliasEquations = []
    defined = Dict()

    nonAliasEquations = []
    findAliases2!(nonAliasEquations, aliases, equations, merge(unknowns, params), defined)

#=
    for eq in equations
        findAliases!(nonAliasEquations, aliases, eq, unknowns, defined)
    end
=#
    loglnModia("\nAlias equations")
    for a in aliases
        loglnModia(prettyfy(a[1]), " = ", prettyfy(a[2]))
    end

    loglnModia("\nAlias definitions")
    for a in aliases
        (nonAlias, state) = getNonAliasVariable(aliases, a[2], unknowns)
#        unknowns[nonAlias.name].state = state  # Temporarily
        loglnModia(prettyfy(a[1]), " := ", prettyfy(nonAlias))
    end
        
    substitutedAliasEquations = []
    for eq in nonAliasEquations
        push!(substitutedAliasEquations, substiteAlias(aliases, eq, unknowns))
    end
    
    equations = substitutedAliasEquations
    flat_model.equations = equations
    
    # Remove alias variables from unknowns
    nonAliased = VariableDict()
    for (n, v) in unknowns
        if !haskey(aliases, GetField(This(), n))
            nonAliased[n] = v
        end
    end
    
    unknowns = nonAliased

    flat_model.variables = copy(unknowns)
    for (p, v) in params
        flat_model.variables[p] = v
    end
    
    if true # PrintFlattened
        loglnModia("\nFlattened model after alias elimination")
        showInstance(flat_model)
        loglnModia()
    end

end

@static if VERSION < v"0.7.0-DEV.2005"
    emptyFieldNames = []
else
    emptyFieldNames = ()
end 

# Type and size deduction

""" 
Modia supports type and size inference, that is, the Variable constructor does not need to specify type and size. However, Pantelides algorithm and removal of singularities require that types and sizes are known.
Types and sizes can be inferred from the start values provided. An outline of such an algorithm is:
- Determine size and type for variables having type, size or start attribute.
- Substitute parameter values in all equations
- Substitute variables with their start values if any in all equations.
- Repeatedly evaluate variables and substitute until no changes
  - If only one remaining unknown, solve for it.
  - Check if the equation only has +/- operators and at least one variable has known size and type. If so, copy size and type to the other variables.
  
Note: This function should be rewritten to ensure consistency, etc.
"""
function deduceVariableAndEquationSizes(flat_model, unknowns, params, equations)
    loglnModia("\nSIZE AND TYPE DEDUCTION")
    # Redefinition of get_start:
    get_start(v::Variable) = v.start
    #  get_start(x) = x

    # Create substitution Dict
    subsValues = Subs() 
    # Create Dicts for variable sizes and types
    varSizes = VariableDict() 
    varTypes = VariableDict() 

    sizeOfType(T) = try
        size(zero(T))
    catch
        nothing
    end

    for (name, var) in unknowns    
        if var.typ != Any && var.size != nothing
            # check consistency!!!
            varSizes[name] = var.size
            varTypes[name] = var.typ
            if extendedLog
                v = name
                loglnModia("    ", v, "[", varSizes[v], "] :: ", varTypes[v], " derived from: T = ", var.typ, " and: size = ", var.size)
            end
        else
            if var.typ != Any
                varTypes[name] = var.typ    
                v = name
                if extendedLog
                    if v in keys(varSizes)
                        loglnModia("    ", v, "[", varSizes[v], "] :: ", varTypes[v], " derived from: T = ", var.typ)
                    else
                        loglnModia("    ", v, "[", "?", "] :: ", varTypes[v], " derived from: T = ", var.typ)
                    end
                end
            end
        
            if var.size != nothing
                varSizes[name] = var.size   
                v = name
                if extendedLog
                    if v in keys(varTypes)
                        loglnModia("    ", v, "[", varSizes[v], "] :: ", varTypes[v], " derived from: size = ", var.size)
                    else
                        loglnModia("    ", v, "[", varSizes[v], "] :: ", "Any", " derived from: size = ", var.size)
                    end
                end            
            end
        end
        # Determine or check size and type for variables having start value
        
        if var.start != nothing
            # check consistency
            siz = size(var.start)
            typ = typeof(var.start)
            if name in keys(varSizes) && siz != varSizes[name]
                println("size($name) and size($name.start) is not consistent.")
            end
        
            #= More elaborate test should be made combining Float constructor with Array{}
            if name in keys(varTypes) && typ != varTypes[name]
                println("typeof($name) and typeof($name.start) is not consistent.")
            end
            =#
            varSizes[name] = size(var.start)
            varTypes[name] = typeof(var.start)
            loglnModia("    ", name, "[", varSizes[name], "] :: ", varTypes[name], " derived from: start = ", var.start)
        end
    end
    #=
    for (k,(name,var)) in enumerate(unknowns)
        s[GetField(This(), name)] = name
    #   s[Der(GetField(This(), name))] = der_name_of(name)
    end
    =#

    # Substitute parameter values
    for (name, var) in params
        subsValues[GetField(This(), name)] = get_value(var)
    end     
    
    subsValues[time_global] = 0.0
    # Substitute variables with their start values if any
    for (name, var) in unknowns
        # st = var.start # get_start(var)
        if var.start != nothing
            subsValues[GetField(This(), name)] = var.start
            subsValues[Der(GetField(This(), name))] = var.start
        end
    end   
    # @show subsValues

    function tryEval(e, eq)
        E = nothing
        try 
            E = eval(e)
        catch err
            if isa(err, ErrorException) && occursin("Unit mismatch", err.msg)
                loglnModia("Warning: ", err.msg, 
                "\n  in expression: ", prettyPrint(e), 
                "\n  in equation:   ", prettyPrint(eq))            
            end
            E
        end
    end

            
    loglnModia("\nEquations")
    equSizes = Dict() # tSize[] # fill(Any(), size(equations))
    equTypes = Dict() # fill(Float64, size(equations))
    equationSubs = copy(equations)
    repeat = true
    
    while repeat
        loglnModia("\nRepeat")
        repeat = false
        for i in 1:length(equations)
            if !haskey(equSizes, i) || !haskey(equTypes, i)
                eq = equations[i]
                # loglnModia()
                # loglnModia(prettyPrint(eq))
                # Substitute parameter values, start values and sofar evaluated variables
                eqsubs = subs(subsValues, eq, false)
                equationSubs[i] = eqsubs

                # Should not use solve due to problem with \ operator finding minimum norm solution for non-square matrix
                vars = []
                findIncidence!(vars, eqsubs)
         
                # If only one remaining unknown, solve for it.
                if false # length(vars) == 1 # disable
                    (e, solved) = SymbolicTransform.solve(eqsubs, vars[1])
                    if solved 
                        # println(prettyPrint(e))
                        eqsubs = e
                    end
                end
                # remove solve
          
                lhs = eqsubs.args[1]
                rhs = eqsubs.args[2]
                
                if typeof(lhs) != GetField && typeof(rhs) == GetField
                    # Swap if equation of form expr = v
                    lhs = eqsubs.args[2]
                    rhs = eqsubs.args[1]
                end
    
                elhs = eq.args[1]
                # investigate equation of type: v = expr or (expr = v due to above swapping)
                if typeof(elhs) == GetField && haskey(varSizes, elhs.name) 
                    if debug
                        @show elhs varSizes[elhs.name]
                    end
                    # Set equation size to the size of the left hand variable
                    
                    equSizes[i] = varSizes[elhs.name]
                    if typeof(rhs) == Expr && rhs.head == :call && rhs.args[1] in [+, -] && length(rhs.args) == 3
                        e1 = rhs.args[2]
                        e2 = rhs.args[3]
                    
                        if typeof(e1) == GetField
                            varSizes[e1.name] = varSizes[elhs.name]
                            v = e1.name
                            if v in keys(varTypes)
                                loglnModia("    ", v, "[", varSizes[v], "] :: ", varTypes[v], " derived from: ", prettyPrint(eq))
                            else
                                loglnModia("    ", v, "[", varSizes[v], "] :: ", "Any", " derived from: ", prettyPrint(eq))                
                            end
                        end
                    
                        if typeof(e2) == GetField
                            varSizes[e2.name] = varSizes[elhs.name]
                            v = e2.name
                            if v in keys(varTypes)
                                loglnModia("    ", v, "[", varSizes[v], "] :: ", varTypes[v], " derived from: ", prettyPrint(eq))
                            else
                                loglnModia("    ", v, "[", varSizes[v], "] :: ", "Any", " derived from: ", prettyPrint(eq))                
                            end
                        end
                    end
                end
          
                if typeof(elhs) == GetField && haskey(varSizes, elhs.name)
                    #=
                    println("\nCHECK OUT")
                    @show prettyPrint(eq)
                    @show elhs.name            
                    @show varSizes[elhs.name]
                    println()
                    =#
                    equSizes[i] = varSizes[elhs.name]
                end
                
                if typeof(elhs) == GetField && haskey(varTypes, elhs.name)
                    #=
                    println("\nCHECK OUT")
                    @show prettyPrint(eq)
                    @show elhs.name            
                    @show varTypes[elhs.name]
                    println()
                    =#
                    equTypes[i] = varTypes[elhs.name]
                end
          
                if typeof(lhs) == GetField 
                    if !haskey(equSizes, i) || !haskey(equTypes, i)
                        RHS = tryEval(rhs, eq)
                        if RHS != nothing && typeof(RHS) != GetField && typeof(RHS) != Der
                            if typeof(RHS) == String || fieldnames(typeof(RHS)) != emptyFieldNames || typeof(RHS) == Tuple{Array{Float64,2},Array{Float64,2}} # Special case for calling qr function
                            # Using isstructtype(typeof(RHS)) did not work as expected so fieldnames(typeof(RHS)) != emptyFieldNames is used instead
                                size_RHS = ()
                            else
                                size_RHS = size(RHS)
                            end
                    
                            loglnModia("  Start values:  ", prettyfy(lhs), " = ", prettyPrint(RHS))
                            subsValues[lhs] = RHS
                            subsValues[Der(lhs)] = RHS
                            varSizes[lhs.name] = size_RHS
                    
                            if !haskey(varTypes, lhs.name) || varTypes[lhs.name] != Any
                                varTypes[lhs.name] = typeof(RHS)
                            end
                    
                            v = lhs.name
                            loglnModia(v, "[", varSizes[v], "] :: ", varTypes[v], " derived from: ", prettyPrint(eq))
                            equSizes[i] = size_RHS
                            equTypes[i] = typeof(RHS)
                            v = lhs.name
                            # loglnModia(v, "[", varSizes[v], "] :: ", varTypes[v])
                            loglnModia("SIZE = ", size_RHS, "\t  TYPE = ", typeof(RHS), "\t   ", prettyPrint(eq))
                            # Since variable value was determined we need to subsitute and repeat checking of equations
                            repeat = true
                        else
#                           loglnModia("Could not evaluate 2: $rhs")
                        end
                    end
                else
                    LHS = tryEval(lhs, eq)
                    RHS = tryEval(rhs, eq)
                    if LHS != nothing && typeof(LHS) != Der && RHS != nothing && typeof(RHS) != Der
                        if typeof(RHS) == String || typeof(RHS) != AbstractArray || fieldnames(typeof(RHS)) != emptyFieldNames
                            size_RHS = ()
                        else
                            size_RHS = size(RHS)
                        end
#=
                        if typeof(LHS) == String || typeof(LHS) != AbstractArray || fieldnames(typeof(LHS)) != emptyFieldNames
                            size_LHS = ()
                        else
                            size_LHS = size(LHS)
                        end
=#
                        size_LHS = size(LHS)
                        if size_LHS != size_RHS 
                            loglnModia("Warning: Not equal size of left and right hand side in equation $(prettyPrint(eq)): $LHS = $RHS")
                        end
                      
                        e = promote(LHS, RHS)
                        t = typeof(e[1])
                        # loglnModia("ASSERT:  ", typeof(LHS), " == ", typeof(RHS))
                        loglnModia("SIZE = ", size_LHS, "  \tTYPE = ", t, "  \t", prettyPrint(eq))
                        equSizes[i] = size_LHS
                        equTypes[i] = t
                        equationSubs[i] = eqsubs
                        #= Fix for call back function
                    elseif LHS != nothing && typeof(LHS) != Der
                        t = typeof(LHS)
                        loglnModia("SIZE = ", size(LHS), "  \tTYPE = ", t, "  \t", prettyPrint(eq))
                        equSizes[i] = size(LHS)
                        equTypes[i] = t
                        equationSubs[i] = eqsubs
                        # else RHS
                        =#
                    end
                end
            end
        end
    end
    
    # Log results
    unknownSizes = false   
    ESizes = [] # tSize[]
    ETypes = []
    loglnModia()
    # @show equSizes
    
    for i in 1:length(equations)
        if !haskey(equSizes, i) || !haskey(equTypes, i)
            eq = equations[i]
            eqs = equationSubs[i]
    
            if !unknownSizes
                loglnModia("Unknown equation sizes or types (scalar Float assumed):")
            end
    
            loglnModia("SIZE = ?  \tTYPE = ?  \t", prettyPrint(eq), ",   \t", prettyPrint(eqs))
            unknownSizes = true
            # Patch waiting for inference in flow equations!!!
            # Assume scalar Float
            push!(ESizes, ()) ### Patch
            push!(ETypes, Float64)  ### Patch
        else
            push!(ESizes, equSizes[i])
            push!(ETypes, equTypes[i])
        end
    end
    # Patch waiting for inference in flow equations!!!
    unknownSizes = false   ### Patch
    
    loglnModia("\nVariables")
    for v in unknowns.keys
        if haskey(varSizes, v) && haskey(varTypes, v)
            loglnModia(v, "[", varSizes[v], "] :: ", varTypes[v])
            
            if flat_model.variables[v].typ == Any
                flat_model.variables[v].typ = varTypes[v]
            end
            
            if flat_model.variables[v].size == nothing 
                flat_model.variables[v].size = varSizes[v]
            end
            
            if flat_model.variables[v].start == nothing
                if varTypes[v] == String
                    z = ""
                    flat_model.variables[v].start = fill(z, varSizes[v])
                else
                    z = 0.0 # zero(varTypes[v]) # Generalize to handle zero for array types
                    flat_model.variables[v].start = if varSizes[v] == (); z else fill(z, varSizes[v]) end
                end
            end

        elseif haskey(varSizes, v)
            loglnModia(v, "[", varSizes[v], "] :: ", "?")
            
            if flat_model.variables[v].size == nothing 
                flat_model.variables[v].size = varSizes[v]
            end
            
            if flat_model.variables[v].start == nothing
                flat_model.variables[v].start = zeros(varSizes[v])  # Generalize
            end

        elseif haskey(varTypes, v)
            if flat_model.variables[v].typ == Any
                flat_model.variables[v].typ = varTypes[v]
            end
            loglnModia(v, "[", "?", "] :: ", varTypes[v])
        
        else
            loglnModia(v, "[?] :: ?")
            unknownSizes = true
        end
    end
    loglnModia()
    
    if unknownSizes
        ModiaLogging.increaseLogCategory(:SizeInferenceError)
        ModiaLogging.closeLogModia()
    end
    @assert !unknownSizes "Not possible to infer size of all variables and equations. Examine log for what variables that need sizes or start values."

    VSizes = [] # tSize[]
    VTypes = []
    
    for v in unknowns.keys
        s = if haskey(varSizes, v); varSizes[v] else () end
        push!(VSizes, s)
        t = if haskey(varTypes, v); varTypes[v] else Any end
        push!(VTypes, t)
    end

    return varSizes, varTypes, equSizes, equTypes, VSizes, VTypes, ESizes, ETypes
end


# -----------------------------------------------------------------
 
function transformStructurally(flat_model)
    equations = flat_model.equations
    params, unknowns = split_variables(vars_of(flat_model))
    setTimeInvariants(params.keys)
  
    loglnModia("Number of parameters: ", length(params))
    loglnModia("Number of unknowns:   ", length(unknowns))
    loglnModia("Number of equations:  ", length(equations))
  
    nFlow = 0
    for (n, v) in unknowns
        if v.flow
            nFlow += 1
        end
    end
    loglnModia("Number of flow variables:  ", nFlow)

    retur = false
    if nFlow == length(equations)
        println("Connector or connector instances")
        ModiaLogging.increaseLogCategory(:ConnectorOrConnectorInstances)
        retur = true
    end

    if flat_model.partial 
        println("Partial model")
        ModiaLogging.increaseLogCategory(:Partial)
        # retur = true
    end

    if retur  
        return nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
    end
  
    #=
    unknownsNames = unknowns.keys
    unknowns_indices = [key => k for (k,key) in enumerate(unknownsNames)] 
    if log 
        printSymbolList("\nunknowns: ", unknownsNames)
    end
    =#

    if aliasElimination
        performAliasElimination!(flat_model, unknowns, params, equations)
    end

    equations = flat_model.equations
    params, unknowns = split_variables(vars_of(flat_model))
  
    loglnModia("Number of parameters: ", length(params))
    loglnModia("Number of unknowns:   ", length(unknowns))
    loglnModia("Number of equations:  ", length(equations))
    
    if deduceSizes
        (varSizes, varTypes, equSizes, equTypes, VSizes, VTypes, ESizes, ETypes) = deduceVariableAndEquationSizes(flat_model, unknowns, params, equations)
    end
  
    #  checkSizes(VSizes, ESizes)
 
    states = Vector()
    deriv = Vector()
    for eq in equations
        findStates!(states, deriv, eq)
    end
    states = unique(states)
    deriv = unique(deriv)
    
    loglnModia("\nINCIDENCE GRAPH")
    unknownsNames = unknowns.keys
    # @show unknownsNames
    unknowns_indices = Dict(key => k for (k, key) in enumerate(unknownsNames))
  
    # Build variable association list. Avar[j] points to entry for derivative of variale j.
    Avar = [findfirst(isequal(GetField(This(), name)), states) for name in unknownsNames]
    #=
    Avar = fill(0, length(unknownsNames))
    for k in 1:length(deriv)
        d = deriv[k]
        j = findfirst(isequal(Symbol(d.base.name)), unknownsNames)  
        Avar[j] = k
    end
    =#
    # Add index offset for deriv vector and append zeros for deriv.
    Avar = [[if a != notFound; a + length(unknownsNames) else 0 end for a in Avar]; fill(0, length(deriv))]

    if false # log 
        printSymbolList("\nUnknowns", unknownsNames, true, true, Avar)
    end

    nonStateVariables = findNonStateVariables(unknowns)
  
    if false
        @show states
        @show deriv
        @show Avar
    end
    #  printSymbolList("\ndifferentiated: ", [states[i].name for i in 1:length(states)])
    printSymbolList("\nNon state variables", nonStateVariables)
  
    statesIndices = []
    realStates = []
    for i in 1:length(states)
        if !(states[i].name in nonStateVariables)
            push!(realStates, states[i])
            push!(statesIndices, unknowns_indices[states[i].name])
        end
    end 

    setRealStates(realStates)
    printSymbolList("Real state variables", realStates)
  
    n = length(unknownsNames) + length(deriv)
    
    equations, equationsIG, variables, assignIG, componentsIG, Avar, Bequ, states, deriv, unassignedNames, incidenceMatrix, VSizes, VTypes, ESizes, ETypes = 
    analyzeStructurally(equations, params, unknowns_indices, deriv, unknownsNames, Avar, statesIndices, states, nonStateVariables, n, realStates, findIncidence!, VSizes, VTypes, ESizes, ETypes, unknowns, flat_model.partial)
  
    if equations == nothing
        return nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing
    end

    #=
    loglnModia("Updated equations")
    for e in equations
        @show e
    end
    =#
  
    return equations, variables, assignIG, componentsIG, Avar, Bequ, states, deriv, unassignedNames, incidenceMatrix, varSizes, varTypes, equSizes, equTypes
end

end