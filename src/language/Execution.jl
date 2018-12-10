"""
Modia module for executing a model including code generation and calling DAE solver.

* Original developer: Toivo Henningsson, Lund
* Developer: Hilding Elmqvist, Mogram AB
* Copyright (c) 2016-2018: Hilding Elmqvist, Toivo Henningsson, Martin Otter
* License: MIT (expat)

"""
module Execution

export simulate_ida
export setOptions

using Base.Meta: quot, isexpr
using DataStructures: OrderedDict

@static if VERSION < v"0.7.0-DEV.2005"
    evaluate(m, x) = eval(m, x)
else
    using LinearAlgebra
    using SparseArrays
    using Dates
    evaluate(m, x) = Core.eval(m, x)
end

import ..Instantiation: Symbolic, Der, Instance, AbstractDict, VariableDict, Variable, Nothing, time_symbol, simulationModel_symbol, vars_of, check_start, GetField, This, time_global, simulationModel_global, eqs_of, get_start, get_dims, model_name_of, operator_table, prettyPrint
import ..Instantiation: Variability, constant, parameter, discrete, continuous

import ModiaMath #0.7
using Unitful
using ..ModiaLogging

const PrintJSONsolved = false 
const showCode = false        # Show the code for initialization and residual calculations
const logComputations = false # Insert logging of variable values in the code
const callF = false
const showJacobian = false

const MODULES = Module[]


global logTiming               # Show timing for each major task, simulate twice to see effect of compilation time.
global storeEliminated
global handleImpulses

function setOptions(options) 
    global storeEliminated = true
    if haskey(options, :storeEliminated)
        global storeEliminated = options[:storeEliminated]
        @show storeEliminated
        delete!(options, :storeEliminated)
    end

    global handleImpulses = false
    if haskey(options, :handleImpulses)
        global handleImpulses = options[:handleImpulses]
        @show handleImpulses
        delete!(options, :handleImpulses)
    end

    global logTiming = false
    if haskey(options, :logTiming)
        global logTiming = options[:logTiming]
        @show logTiming
        delete!(options, :logTiming)
    end
end

# ----------------------------------

function split_variables(src::VariableDict)
    params, vars = VariableDict(), VariableDict()
    for (name, var) in src
        if !isa(var, Instance) && (!isa(var, Variable) || var.variability <= parameter)
            params[name] = var
        else                                                   
            vars[name] = var
        end
    end
    (params, vars)
end

# Substitution of variable names for vector elements, e.g. z => x[3]
const Subs = Dict{Symbolic,Any}

subs(s::Subs, ex, complete::Bool) = get(s, ex, ex)
function subs(s::Subs, ex::Symbolic, complete::Bool)
        if haskey(s, ex);  s[ex]
            elseif complete;   ModiaLogging.increaseLogCategory(:NoSubstitution); error("No substitution provided for: ", ex)
        else               ex
    end
end

function subs(s::Subs, ex::Expr, complete::Bool)
    if isexpr(ex, :quote)
        ex
    else
        Expr(ex.head, [subs(s, arg, complete) for arg in ex.args]...)
    end
end

subs(s::Subs, ex) = subs(s, ex, false)

get_value(v::Variable) = v.value
get_value(x) = x

code_reshape(ex, dims::Int...) = :( $(quot(reshape))($ex, $dims) )

code_state_read(x::Symbol, k::Int) = :($x[$k])
code_state_read(x::Symbol, k::Int, size::Int) = :($x[$k:$(k + size - 1)])
code_state_read(x::Symbol, k::Int, dims::Int...) = code_reshape(code_state_read(x, k, prod(dims)), dims...)

get_name(g::GetField) = (@assert g.base == This(); g.name)

mark_solved_eq(ex, solved::Bool) = error("Only support = and := equations, got ", ex)

function mark_solved_eq(ex::Expr, solved::Bool)
        if isexpr(ex, [:(=), :(:=)], 2);  Expr(solved ? :(:=) : :(=), ex.args...)
        else                              error("Only support = and := equations, got ", ex)
    end
end

mark_solved_equations(eqs::Vector, solved::Vector{Bool}) = [mark_solved_eq(eq, sol) for (eq, sol) in zip(eqs, solved)]

find_diffstates!(d::Set{Symbolic}, ex) = nothing
find_diffstates!(d::Set{Symbolic}, der::Der) = (push!(d, der.base); nothing)
function find_diffstates!(d::Set{Symbolic}, ex::Expr)
    if !isexpr(ex, :quote)
        for arg in ex.args
            find_diffstates!(d, arg)
        end
    end
    nothing
end

function code_eliminated_func(fname, unpack, eliminated_computations, vars, x::Symbol, der_x::Symbol)
    @gensym xs der_xs results ts k
    fname = Symbol(fname)

    alloc_eliminated = []
    push_eliminated = []

    for computation in eliminated_computations
        name = (computation::Expr).args[1]::Symbol
        res_name = gensym(string(name, "_results"))

        var = vars[name]
        T, dims = var.typ, get_dims(var)
        # if T == A; T = Float64; end
        # if T != Any || T != Float64; T = Any; end  # Hack to handle units
        if typeof(T) <: Unitful.Unitlike; T = Float64; end # To handle units
        # @show T, string(var)
        
        if !isempty(dims);  T = Array{T,length(dims)};  end

        @static if VERSION < v"0.7.0-DEV.2005"
            push!(alloc_eliminated, :($res_name = $results[$(string(name))] = Vector{$(quot(T))}(0)))
        else
            push!(alloc_eliminated, :($res_name = $results[$(string(name))] = Vector{$(quot(T))}(undef, 0)))
        end
        push!(push_eliminated,  :($(quot(push!))($res_name, $name)))
    end
    # @show eliminated_computations
    # @show push_eliminated
    
    eliminated_code = quote
        function $(fname)($results, $ts, $xs, $der_xs)
            $(alloc_eliminated...)
            for ($k, $time_symbol) in enumerate($ts)
#                $x = $(quot(vec))($xs[$k, :])
                # For 0.7:
                $x = $xs[$k, :]
#                $der_x = $(quot(vec))($der_xs[$k, :])
                $der_x = $der_xs[$k, :]
                $(unpack...)
                $(eliminated_computations...)
                $(push_eliminated...)
            end
            nothing
        end
    end

    if showCode
        @show eliminated_code
    end

    eliminated_code
end

der_name_of(name::Symbol) = Symbol("der(", name, ")")

function store_eliminated!(results::Vector{Vector}, k::Int, result)
    dest = results[k]
    if !isa(result, eltype(dest))
        T = promote_type(eltype(dest), typeof(result))
        dest = results[k] = Vector{T}(dest)
    end
    push!(dest, result)
end

function substituteExpr(ex, equations, s)
    if true # Check expression
        for eq in equations
            if isexpr(eq, :(:=), 2)
                # lhs_name, rhs = get_name(eq.args[1]), subs(s, eq.args[2], true)
                if ex == eq.args[1]
                    ex = eq.args[2]
                    break
                end
            end            
        end
    end
    cond = subs(s, ex, true)
    evalCond = eval(cond)
    return evalCond
end

global oldCond = false
global F_Dict = Dict()

struct Eval{F}
    f::F
end
(m::Eval)(x...) = eval(m.f)(x...)


@static if VERSION < v"0.7.0-DEV.2005"
    const letArgs = 1
else
    const letArgs = 2
end      
      
function prepare_ida(instance::Instance, first_F_args, initial_bindings::AbstractDict{Symbol,Any}; store_eliminated=false, need_eliminated_f=false)
    proceed = zeros(Bool, 1)
    global F_Dict

    params, vars = split_variables(vars_of(instance))
    for (name, var) in vars
        check_start(var, name)
    end

    s = Subs()
    for (k, (name, var)) in enumerate(vars)
        s[GetField(This(), name)] = name
        s[Der(GetField(This(), name))] = der_name_of(name)
    end

    for (name, var) in params
        s[GetField(This(), name)] = get_value(var)
    end

    s[time_global] = time_symbol
    s[simulationModel_global] = simulationModel_symbol

    computations = []
    # push!(computations, :(if ModiaMath.ModiaToModiaMath.isPreInitial($simulationModel_symbol)
    # $simulationModel_symbol.store = Synchronous.Store()
    # end))
    if logComputations
        push!(computations, :(if !$proceed[1]; println("\nCOMPUTATIONS") end))
        # push!(computations, :(@show $simulationModel_symbol))
    end
    eliminated_computations = []

    residuals = Symbol[]
    eliminated = Symbol[]
    states = copy(vars)

    function getType(vars, v)
        vars[v].typ
    end

    modeConditions = []
    for eq in eqs_of(instance)
        if isexpr(eq, :(:=), 2)
            lhs_name, rhs = get_name(eq.args[1]), subs(s, eq.args[2], true)

            @assert haskey(states, lhs_name)
            delete!(states, lhs_name)
            T = getType(vars, lhs_name)
            ex = :($lhs_name = $rhs)
            # ex = :($lhs_name::$T = $rhs)
            push!(computations, ex)
            if logComputations
                push!(computations, :(if !$proceed[1]; @show typeof($lhs_name) end))
                push!(computations, :(if !$proceed[1]; @show $lhs_name end))
            end
            push!(eliminated_computations, ex)
            push!(eliminated, lhs_name)
        elseif isexpr(eq, :(=), 2)
            eq = subs(s, eq, true)
            lhs, rhs = eq.args
            residual = gensym("residual")
            push!(residuals, residual)
            push!(computations, :($residual = $rhs - $lhs))
            # push!(computations, :($residual = Unitful.ustrip($rhs - $lhs)) )  # need using Unitful in function
            if logComputations
                push!(computations, :(if !$proceed[1]; @show $residual end))
            end

        else # conditional equations
            time = initial_bindings[time_symbol]
            s[time_global] = time
            condExpr = eq.args[1]
            evalCond = substituteExpr(condExpr, eqs_of(instance), s)
            push!(modeConditions, evalCond)
            if time == 0.0
                println("Mode condition: ", prettyPrint(condExpr), " is $evalCond at time=$time")
                oldCond = evalCond
            end

            global oldCond
            if evalCond != oldCond
                println("Mode change: ", prettyPrint(condExpr), " became ", evalCond, " at time ", time)
            end
            oldCond = evalCond

            e = if evalCond; eq.args[2] else eq.args[3] end
            e = subs(s, e, true)

            eblock = e.args

            for ei in eblock
                if ei != nothing
                    lhs, rhs = ei.args
                    residual = gensym("residual")
                    push!(residuals, residual)
                    push!(computations, :($residual = $rhs - $lhs))
                    if logComputations
                        push!(computations, :(if !$proceed[1]; @show $residual end))
                    end
                end
            end
                        
            # error("Unsupported equation type: ", eq)
        end
    end

    d = Set{Symbolic}()
    for eq in eqs_of(instance);  find_diffstates!(d, eq);  end

    x0 = Float64[]
    x_nominal = Float64[]
    diffstates = Bool[]
    for (name, var) in states
        is_diffstate = GetField(This(), name) in d
        s = get_start(var)
        if isa(s, AbstractArray)
            append!(x0, vec(s))
            append!(x_nominal, vec(var.nominal === nothing ? fill(1.0, var.size) : var.nominal))
            # @show name, vec(s)            
            append!(diffstates, fill(is_diffstate, length(s)))
        else
            push!(x0, s)
            push!(x_nominal, var.nominal === nothing ? 1.0 : var.nominal)
            # push!(x0, ustrip(s))
            # @show name, s  
            push!(diffstates, is_diffstate)
        end
    end
    # @show x0 size(x0)
    # @show x_nominal
    # @show diffstates
    #=
    println("State vector allocation:")
    i = 0
    for ((name, var), offset) in zip(states, state_offsets)
        der_name = der_name_of(name)
        print("$name = x[$offset]")
        i += 1
        if diffstates[i]
          println("   $der_name = der_x[$offset]")
        else
          println()
        end
    end
    =#    

    # Create mapping between states and state vector
    state_sizes = [prod(get_dims(var)) for var in values(states)]
    state_offsets = cumsum(vcat([1], state_sizes[1:end - 1]))

    if state_sizes != Any[]
        state_size = sum(state_sizes)
    else
        state_size = 0
    end

    if state_size > 0
        ModiaLogging.increaseLogCategory(:DynamicModel)
    else
        ModiaLogging.increaseLogCategory(:StaticModel)
    end

    loglnModia("statesize = ", state_size)
    
    # if ! haskey(F_Dict, modeConditions)

    @gensym x der_x r
    unpack = []

    if logComputations
        push!(unpack, :(
          if !$proceed[1]; 
            println("Press enter to continue, q to stop, p to proceed: "); l = readline(stdin); if l != "" && l[1] == 'q'; error("quit") elseif l != "" && l[1] == 'p'; $proceed[1] = true end
        end))
        push!(unpack, :(if !$proceed[1]; println("\nUnpack:") end))
    end
      
    # println("State vector allocation:")
    initials = [] 
    i = 1
    for ((name, var), offset) in zip(states, state_offsets)
        der_name = der_name_of(name)
          
        if true # diffstates[i]  # Wonder about reason???
            push!(unpack, :(global $name = $(code_state_read(x, offset, get_dims(var)...))))
            push!(unpack, :($der_name = $(code_state_read(der_x, offset, get_dims(var)...))))
        else
            push!(unpack, :(global $name = $(code_state_read(der_x, offset, get_dims(var)...))))
        end      
        i += prod(get_dims(var))        
        if logComputations
            push!(unpack, :(if !$proceed[1]; @show typeof($name), typeof($der_name) end))        
            push!(unpack, :(if !$proceed[1]; @show $name, $der_name end))
        end

        push!(initials, :($name = $(quot(get_start(var)))))
        if logComputations
            push!(initials, :(@show typeof($name) $name) )
        end
        # push!(initials, :($der_name = 0 * $(quot(get_start(var))) / SIUnits.Second) )
        push!(initials, :($der_name = 0 * $(quot(get_start(var))) ))
    end

    for (name, value) in initial_bindings
        push!(initials, :($name = $(quot(value))))
    end

    initial_ex = :(let; end)
    initial_body = (initial_ex.args[letArgs].args)::Vector{Any}
    @assert isempty(initial_body)

    append!(initial_body, initials)
    append!(initial_body, instance.initial_pre)
    append!(initial_body, computations)
    append!(initial_body, instance.initial_post)
    # todo: avoid ... here too?
    push!(initial_body, :(Any[$(residuals...)], Any[$(eliminated...)]))

    if showCode
        println("\nINITIALIZATION CODE")
        @show initial_ex
    end
    
    initial_residuals, initial_eliminated = eval(initial_ex)
    residual_dims = [size(res) for res in initial_residuals]
    # @show residual_dims
    residual_sizes = [prod(dims) for dims in residual_dims]
    # @show residual_sizes
    
    if residual_sizes != Any[]
        residual_size = sum(residual_sizes)
    else
        residual_size = 0
    end
    
    loglnModia("residual_size = ", residual_size)
    @assert state_size == residual_size
      
    residual_offsets = cumsum(vcat([1], residual_sizes[1:end - 1]))
    # @show residual_offsets

    pack = []
    for (residual, dims, offset) in zip(residuals, residual_dims, residual_offsets)
        rhs = residual
        if dims === ()
            lhs = code_state_read(r, offset)
        else
            sz = prod(dims)
            lhs = code_state_read(r, offset, sz)
            if length(dims) > 1;  rhs = code_reshape(rhs, sz);  end
        end
        push!(pack, :($lhs = $rhs))
    end

    store_elim = :nothing
    if store_eliminated
        @gensym eliminated_results
        store_elims = [
              :( $(quot(store_eliminated!))($eliminated_results, $k, $name) )
              for (k, name) in enumerate(eliminated)
          ]
        store_elim = quote
            if ($(quot(is_output_point))($simulationModel_symbol))
                $eliminated_results = $(quot(get_eliminated_results_object))($simulationModel_symbol)
                $(store_elims...)
            end
        end
    end

    name = Symbol("F_", model_name_of(instance))
    F_code = :( function $(name)($(first_F_args...), _t, $x, $der_x, $r, _w)
#                  $r[0] = 0.0 # For the case of no differential equations
                  end )
    F_body = (F_code.args[2].args)::Vector{Any}
    # @assert isempty(F_body)

    append!(F_body, unpack)
    push!(F_body, :($time_symbol = $simulationModel_symbol.simulationState.time)) ####might break Sundials (not DAE)
    append!(F_body, instance.F_pre)
    append!(F_body, computations)
    append!(F_body, instance.F_post)
    append!(F_body, pack)
    # push!(  F_body, store_elim)
    push!(F_body, nothing)
    F_code = quote
        let
            $(F_code)
        end
    end
      
    if showCode
        println("\nFUNCTION F CODE")
        @show F_code
    end

    # F = Eval(F_code)
    push!(MODULES, Module())
    index_Module_F = length(MODULES)
    F = evaluate(MODULES[index_Module_F], F_code)
    #=
      if modeConditions != []
        F_Dict[modeConditions] = (F, initial_eliminated)
      end
    else
      (F, initial_eliminated) = F_Dict[modeConditions]
    end
    =#   

    if need_eliminated_f
        eliminated_code = code_eliminated_func(string("eliminated_", model_name_of(instance), "!"),
            unpack, eliminated_computations, vars, x, der_x)
        eliminated_f = evaluate(MODULES[index_Module_F], eliminated_code)
    else
        eliminated_f = nothing
    end

    eliminated_Ts = OrderedDict{Symbol,Type}()
    for (name, value) in zip(eliminated, initial_eliminated)
        eliminated_Ts[name] = typeof(value)
    end

    der_x0 = zeros(size(x0)) # todo: allow other initial guesses for der_x?
    (F, eliminated_f, x0, der_x0, x_nominal, diffstates, params, states, state_sizes, state_offsets, eliminated_Ts)
end

# Temporary until we can store it in SimulationModel
global_elim_results = Vector{Vector}()

# Temporary definitions
import ModiaMath.ModiaToModiaMath.ModiaSimulationModel
#is_output_point(m::ModiaSimulationModel) = true # Should not always be true!
is_output_point(m::ModiaSimulationModel) = false # Should not always be true!
get_eliminated_results_object(m::ModiaSimulationModel) = global_elim_results

function extract_results_ida(x_res, der_x_res, states, state_offsets, params)
    results = Dict{AbstractString,Any}()
    # @show keys(states)
    for ((name, var), offset) in zip(states, state_offsets)
        inds = get_dims(var) == () ? offset : (offset:(offset + prod(get_dims(var)) - 1))
        results[string(name)] = x_res[:,inds]
        results[string("der(", name, ")")] = der_x_res[:,inds]
    end

    for (name, par) in params
        results[string(name)] = get_value(par)
    end
    results
end


# -----------------------------

function simulate_ida(instance::Instance, t, args...; kwargs...)
    simulate_ida(instance, collect(Float64, t), args...; kwargs...)
end

function simulate_ida(instance::Instance, t::Vector{Float64},
                      jac::Union{SparseMatrixCSC{Bool,Int},Nothing};# =nothing;
                      log=false, relTol=1e-4, hev=1e-8, maxSparsity=0.1,
                      store_eliminated=storeEliminated)

    if PrintJSONsolved
        printJSONforSolvedEquations(instance)
    end
    
    initial_bindings = Dict{Symbol,Any}(time_symbol => t[1])
    initial_m = ModiaSimulationModel()

    initial_bindings[simulationModel_symbol] = initial_m
    
    prep = prepare_ida(instance, [simulationModel_symbol], initial_bindings, store_eliminated=storeEliminated, need_eliminated_f=storeEliminated)
    F, eliminated_f, x0, der_x0, x_nominal, diffstates, params, states, state_sizes, state_offsets, eliminated = prep
    
    eliminated_results = Vector[Vector{T}() for T in values(eliminated)]
    # temporary until we can store it in simulationModel
    global global_elim_results = eliminated_results

    if callF || handleImpulses || showJacobian
        callResidualFunction(F, callF, handleImpulses, showJacobian, x0, der_x0, diffstates, instance)
    end
    
    xNames = fill("[]", length(x0))    
    ii = 0

    for (name, var) in states
        ii += 1
        name = string(name)
        xNames[state_offsets[ii]] = name

        if state_sizes[ii] > 1
            dims = get_dims(var)
            dimsArray = collect(dims)

            if length(dimsArray) == 1
                for j in 1:dimsArray[1]
                    xNames[state_offsets[ii] + j - 1] = name * "[" * string(j) * "]"          
                end
            else
                xNames[state_offsets[ii]] = name * "[]"
                xNames[state_offsets[ii] + state_sizes[ii] - 1] = string(collect(dimsArray))
            end
        end
    end
    # @show xNames
    getVariableName(model::Any, vcat::ModiaMath.DAE.VariableCategory, vindex::Int) = vcat == ModiaMath.DAE.Category_X ? xNames[vindex] :
                                                                                    (vcat == ModiaMath.DAE.Category_W ? "w[" * string(vindex) * "]" :
                                                                                                                        "der(" * xNames[vindex] * ")" )

    start = now()
    
    if length(x0) > 0    
        if false
            m = ModiaSimulationModel(model_name_of(instance), F, x0, der_x0, jac;
                        xw_states=diffstates, maxSparsity=maxSparsity, nc=1, hev=hev, nz=initial_m.nz_preInitial,
                        xNames=xNames, x_nominal=x_nominal)
        else
            m = ModiaSimulationModel(string(model_name_of(instance)), F, x0, getVariableName; x_nominal=x_nominal,
                        maxSparsity=maxSparsity, nc=1, nz=initial_m.nz_preInitial, hev=hev, jac=jac, x_fixed=diffstates)
        end 

        if logTiming
            print("\n  ModiaMath:           ")
            @time ModiaMath.ModiaToModiaMath.simulate(m, t; log=log, tolRel=relTol)
            # @show now()-start
            if logTiming
                print("  ModiaMath:           ")
                @time ModiaMath.ModiaToModiaMath.simulate(m, t; log=log, tolRel=relTol)
                # @show now()-start
            end
        else
            ModiaMath.ModiaToModiaMath.simulate(m, t; log=log, tolRel=relTol)
        end

        (t_res, x_res, der_x_res) = ModiaMath.ModiaToModiaMath.getRawResult(m)
        # @show now()-start
    else
        t_res = t
        x_res = t
        der_x_res = t
    end

    results = extract_results_ida(x_res, der_x_res, states, state_offsets, params)
    # @show keys(results)
    # @show now()-start
    if store_eliminated
        Base.invokelatest(eliminated_f, results, t_res, x_res, der_x_res)
        # @show now()-start
    else
        for (name, result) in zip(keys(eliminated), eliminated_results)
            results[string(name)] = result
        end
        # @show now()-start
    end
    results["time"] = t_res
    results
end

# -----------------------------


"""
Experimental function for calculating rank of system matrix and for impulse handling
"""
function callResidualFunction(F, callF, handleImpulses, showJacobian, x0, der_x0, diffstates, instance)
    n = length(x0)
    mF = ModiaSimulationModel()
    r0 = fill(0.0, n)
    w = []

    Base.invokelatest(F, mF, 0, x0, der_x0, r0, w)
    # @show x0 der_x0 r0

    r = fill(0.0, n)
    delta = 1E-10
    A = fill(0.0, (n, n))
    independent = 0

    for i in 1:n
        X0 = copy(x0)
        X0[i] += delta
        Base.invokelatest(F, mF, 0, X0, der_x0, r, w)
        delta_r = r - r0
        A[:,i] = delta_r / delta

        if rank(A[:, 1:i]) < independent + 1
            if showJacobian
                println("Dependent column: ", i)
                if all(delta_r .== 0)
                    println("Zero disturbance in residuals")
                else
                    @show delta_r
                    x = A[:, 1:i - 1] \ delta_r
                    y = [if abs(xx) < 0.1; 0.0 else xx end for xx in x]
                    nonzeros = [i for i in 1:length(y) if y[i] != 0.0]
                    # @show x y nonzeros
                end
            end
        else
            independent += 1      
        end
    end

    der_index = [i for i in 1:length(diffstates) if diffstates[i]]
    E = fill(0.0, (n, n))
    independent = 0

    for i in 1:n
        if i in der_index
            der_X0 = copy(der_x0)
            der_X0[i] += delta
            Base.invokelatest(F, mF, 0, x0, der_X0, r, w)
            delta_r = r - r0
            E[:,i] = delta_r / delta

            if rank(E[:, 1:i]) < independent + 1
                if showJacobian
                    println("Dependent column: ", i)
                    if all(delta_r .== 0)
                        println("Zero disturbance in residuals")
                    end
                end
            else
                independent += 1      
            end
        end
    end

    if showJacobian
        @show A
        @show size(A)
        @show rank(A)
        @show der_index
        @show E
        @show size(E)
        @show rank(E)
        J = A + E / delta
        @show rank(J)
        if rank(J) != size(A, 1)
            println("The implicit system matrix does not have full rank.")
        end
    end

    if handleImpulses
        oldx0 = copy(x0)
        X0 = copy(x0)
        der_X0 = copy(der_x0)
        repeat = true
        while repeat
            # Make Implicit Euler step to initialize x0
            J = A + E / delta
            # J*dx + r0 = 0
            dx = J \ -r0
            X0 += dx
            der_X0 = (X0 - oldx0) / delta
            Base.invokelatest(F, mF, 0, X0, der_X0, r0, w)
            repeat = norm(dx) > 1E-3  # Change to r0 after updating der_x0?
        end

        params, vars = split_variables(vars_of(instance))
        names = collect(keys(vars))
        #=    
        for i in 1:length(x0)
            if abs(X0[i]-oldx0[i]) > 1E5
                println("Dirac impulse in variable ", names[i])
            end        
        end
        =#
        for i in der_index
            if abs(X0[i] - oldx0[i]) > 1E-7
                println("Dirac impulse caused a change in the state ", names[i], " from $(oldx0[i]) to $(X0[i])")
            end
        end

        oldx0 = copy(x0)
        repeat = true
        Base.invokelatest(F, mF, 0, x0, der_x0, r0, w)
        while repeat
            # Make Implicit Euler step to initialize x0
            J = A + E / delta
            # J*dx + r0 = 0
            dx = J \ -r0
            x0 += dx
            Base.invokelatest(F, mF, 0, x0, der_x0, r0, w)
            repeat = norm(dx) > 1E-3  # Change to r0 after updating der_x0?
        end
    end
end

# -------------------------

#= causes warning
import JSON
JSON.lower(t::This) = "this."
JSON.lower(s::SIUnits.SIUnit) = string(s)
=#

"""
Experiemental code for printing the AST of solved equations in JSON format
"""
function printJSONforSolvedEquations(instance)
    modelName = string(model_name_of(instance))
    fileName = Base.Filesystem.realpath(modelName * ".json")
    println("\nWriting model in JSON format to file \"", fileName, "\"")
    indent = 3

    open(fileName, "w") do file
        variables = JSON.json(instance.variables, indent)
        print(file, "{\n")
        print(file, "\"name\":\"", modelName, "\",\n")
        print(file, "\"variables\":\n" * variables)

        jsonSolved = JSON.json(makeJSON(instance.equations), indent)
        jsonSolved = replace(jsonSolved, "\\", "")
        print(file, ",\n")
        print(file, "\"equations\":\n", jsonSolved)
        print(file, "}\n")
    end

    # jsonSolved = JSON.json(makeJSON(instance))
    # jsonSolved = replace(jsonSolved, "\\", "")
    # println(jsonSolved)
end



#makeJSON(der::Der) = Symbol("der("*string(der.base.name)*")")
makeJSON(der::Der) = Symbol("der(" * replace(string(der.base.name), ".", "_") * ")")
#makeJSON(get::GetField) = get.name
makeJSON(get::GetField) = Symbol(replace(string(get.name), ".", "_")) # get.name # Handle dummy derivatives
makeJSON(s::Symbol) = string(s)

#makeJSON(ex) = ex
makeJSON(ex) = get(operator_table, string(ex), string(ex))
#makeJSON(s::Instance) = s
function makeJSON(ex::Array{Any})
    [makeJSON(e) for e in ex]
end

function makeJSON(ex::Expr)
    if isexpr(ex, :quote) || isexpr(ex, :line)
        nothing
    elseif isexpr(ex, :block)
        makeJSON(ex.args[2])
    else
        Expr(ex.head, [makeJSON(arg) for arg in ex.args]...)
    end
end

end
