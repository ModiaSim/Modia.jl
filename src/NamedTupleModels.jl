#=
Handles models defined as named tuples.

* Developer: Hilding Elmqvist, Mogram AB
* First version: January 2021
* License: MIT (expat)

=#
export mergeModels, recursiveMerge, Redeclare, showModel, @showModel, Model, Map, Par, Var, setLogMerge,
    constant, parameter, input, output, potential, flow, interval, @info_str

using Base.Meta: isexpr
using DataStructures: OrderedDict
using Unitful
using ModiaBase.Symbolic: removeBlock, prepend

function stringifyDefinition(v)
    if typeof(v) in [Symbol, Expr]
        v = removeBlock(v)
        v = ":(" * string(v) * ")"
    end
    return v
end

function showModel(m, level=0)
    level += 1
#    print("  "^(level-1))
    if typeof(m) <: NamedTuple && haskey(m, :class)
        print(m.class)
    end
    println("(")
    for (k, v) in zip(keys(m), m)
        if typeof(v) <: NamedTuple
            print("  "^level, k, " = ")
            showModel(v, level)
        elseif k != :class
            println("  "^level, k, " = ", stringifyDefinition(v), ",")
        end
    end
    println("  "^(level-1), "),")
end

macro showModel(model)
    modelName = string(model)
    mod = :( print($modelName, " = "); showModel($model); println() )
    return esc(mod)
end

global logMerge = false

function setLogMerge(val)
    global logMerge
    logMerge = val
end

function mergeModels(m1::NamedTuple, m2::NamedTuple, env=Symbol())
    mergedModels = OrderedDict{Symbol,Any}(pairs(m1))
    for (k,v) in collect(pairs(m2))
        if typeof(v) <: NamedTuple
            if k in keys(mergedModels) && ! (:_redeclare in keys(v))
                if logMerge; print("In $k: ") end
                m = mergeModels(mergedModels[k], v, k)
                mergedModels[k] = m
            elseif :_redeclare in keys(v)
                if logMerge; println("Redeclaring: $k = $v") end
                mergedModels[k] = v
            elseif nothing in values(v) # TODO: Refine

            else
                if !(:_redeclare in keys(mergedModels))
                    if logMerge; print("Adding: $k = "); showModel(v, 2) end
                end
                mergedModels[k] = v
            end
        elseif v === nothing
            if logMerge; println("Deleting: $k") end
            delete!(mergedModels, k)
        else
            if logMerge
                if k in keys(mergedModels)
                    if mergedModels[k] != v
                        println("Changing: $k = $(stringifyDefinition(mergedModels[k])) to $k = $(stringifyDefinition(v))")
                    end
                elseif !(:_redeclare in keys(mergedModels))
                    println("Adding: $k = $(stringifyDefinition(v))")
                end
            end
            mergedModels[k] = v
        end
    end
#    delete!(mergedModels, :class)
    return (; mergedModels...) # Transform OrderedDict to named tuple
end

Model(; kwargs...) = (; class = :Model, kwargs...)

Map(; kwargs...) = (; class = :Map, kwargs...)

Par(; kwargs...) = Map(; class = :Par, kwargs...)

Var(;kwargs...) = (; class = :Var, kwargs...)
Var(value; kwargs...) = (; class = :Var, value=value, kwargs...)

Redeclare = ( _redeclare = true, )

constant = Var(constant = true)
parameter = Var(parameter = true)
input = Var(input = true)
output = Var(output = true)
potential = Var(potential = true)
flow = Var(flow = true)

interval(min, max) = Var(min=min, max=max)

macro info_str(text)
    Var(info=text)
end

Base.:|(m::NamedTuple, n::NamedTuple) =  mergeModels(m, n) # TODO: Chane to updated recursiveMerge
Base.:|(m, n) = if !(typeof(n) <: NamedTuple); recursiveMerge(m, (; value=n)) else recursiveMerge(n, (value=m,)) end

recursiveMerge(x, ::Nothing) = x
recursiveMerge(x, y) = y
#recursiveMerge(x::Expr, y::Expr) = begin dump(x); dump(y); Expr(x.head, x.args..., y.args...) end
recursiveMerge(x::Expr, y::Tuple) = begin x = copy(x); xargs = x.args; xargs[y[2]] = y[3]; Expr(x.head, xargs...) end

function recursiveMerge(nt1::NamedTuple, nt2::NamedTuple)
    all_keys = union(keys(nt1), keys(nt2))
#    all_keys = setdiff(all_keys, [:class])
    gen = Base.Generator(all_keys) do key
        v1 = get(nt1, key, nothing)
        v2 = get(nt2, key, nothing)
        key => recursiveMerge(v1, v2)
    end
    return (; gen...)
end


function unpackPath(path, sequence)
    if typeof(path) == Symbol
        push!(sequence, path)
    elseif isexpr(path, :.)
        unpackPath(path.args[1], sequence)
        push!(sequence, path.args[2].value)
    elseif isexpr(path, :ref)
        unpackPath(path.args[1], sequence)
        push!(sequence, path.args[2:end]...)
    end
end

function collectConnector(model)
    potentials = []
    flows = []
    for (k,v) in zip(keys(model), model)
        if typeof(v) <: NamedTuple && :class in keys(v) && v.class == :Var ||
            typeof(v) <: NamedTuple && :variable in keys(v) && v.variable
            if :potential in keys(v) && v.potential
                push!(potentials, k)
            end
            if :flow in keys(v) && v.flow
                push!(flows, k)
            end
        end
    end
    return potentials, flows 
end

function convertConnections!(connections, model, modelName, logging=false)
#    println("\nconvertConnections")
#    showModel(model)
#    @show connections
#    println()
    connectEquations = []
    alreadyConnected = []
    for i in 1:length(connections)
        c = connections[i]

        if c.head == :tuple
            connected = c.args

            potential1 = nothing
            inflow = 0
            outflow = 0
            signalFlow1 = nothing
            connectedOutput = nothing
            potentials1 = []
            flows1 = []
            for con in connected
                if con in alreadyConnected
                    error("Already connected: $con, found in connection: $connected")
                end
                push!(alreadyConnected, con)

                sequence = []
                unpackPath(con, sequence)
                mod = model
                for s in sequence[1:end-1]
                    if s in keys(mod)
                        mod = mod[s]
                    else
                        error("Invalid path $con: $s not found in $(keys(mod))")
                    end
                end

                if sequence[end] in keys(mod)
                    mod = mod[sequence[end]]

                    if :input in keys(mod) && mod.input || :output in keys(mod) && mod.output
                        signalFlow = con
                        if signalFlow1 !== nothing
                            push!(connectEquations, :($signalFlow1 = $signalFlow))
                        end
                        signalFlow1 = signalFlow

                        if :output in keys(mod) && mod.output
                            if connectedOutput != nothing # TODO: refine logic concerning inner and outer inputs and outputs
#                                error("It is not allowed to connect two outputs: ", connectedOutput, ", ", con)
                            end
                            connectedOutput = con
                        end
                    else
                        potentials, flows = collectConnector(mod)
                        # Deprecated
                        if :potentials in keys(mod)
                            potentials = vcat(potentials, mod.potentials.args)
                        end
                        if :flows in keys(mod)
                            flows = vcat(flows, mod.flows.args)
                        end

                        if length(potentials1) > 0 && potentials != potentials1
                            error("Not compatible potential variables: $potentials1 != $potentials, found in connection: $connected")
                        end
                        for p in potentials
                            potential = append(con, p)
                            if potential1 != nothing
                                push!(connectEquations, :($potential1 = $potential))
                            end
                            potential1 = potential
                        end
                        potentials1 = potentials

                        if length(flows1) > 0 && flows != flows1
                            error("Not compatible flow variables: $flows1 != $flows, found in connection: $connected")
                        end
                        for f in flows
                            flowVar = append(con, f)
                            if length(sequence) == 1
                                if inflow == 0
                                    inflow = flowVar
                                else
                                    inflow = :($inflow + $flowVar)
                                end
                            else
                                if outflow == 0
                                    outflow = flowVar
                                else
                                    outflow = :($outflow + $flowVar)
                                end
                            end
                        end
                        flows1 = flows
                    end
                else
                    # Deprecated
                    signalFlow = con
                    if signalFlow1 !== nothing
                        push!(connectEquations, :($signalFlow1 = $signalFlow))
                    end
                    signalFlow1 = signalFlow
                end
            end
            if inflow != 0 || outflow != 0
                push!(connectEquations, :($inflow = $outflow))
            end
        end
    end
    if length(connectEquations) > 0 && logging
        println("Connect equations in $modelName:")
        for e in connectEquations
            println("  ", e)
        end
    end
    return connectEquations
end

function flattenModelTuple!(model, modelStructure, modelName; unitless = false, log=false)
    connections = []
    extendedModel = merge(model, NamedTuple())
    for (k,v) in zip(keys(model), model)
#        @show k v typeof(v)
        # Deprecated
        if k in [:inputs, :outputs, :potentials, :flows, :class]
            if k in [:inputs, :outputs, :potentials, :flows]
                printstyled("  Deprecated construct in model $modelName: ", color=:yellow)
                println("$k = $v\n    Use: ... = $(string(k)[1:end-1]) ...")
            end
        elseif k == :init
            printstyled("  Deprecated construct in model $modelName: ", color=:yellow)
            println("$k = $v\n    Use: ... = Var(init=...) ...")
            for (x,x0) in zip(keys(v), v)
                if x != :class
                    if unitless && typeof(x0) != Expr
                        x0 = ustrip(x0)
                    end
                    modelStructure.init[x] = x0
                    modelStructure.mappedParameters = (;modelStructure.mappedParameters..., x => x0)
                end
            end
        elseif k == :start
            printstyled("  Deprecated construct in model $modelName: ", color=:yellow)
            println("$k = $v\n    Use: ... = Var(start=...) ...")
            for (s,s0) in zip(keys(v), v)
                if s != :class
                    if unitless
                        s0 = ustrip(s0)
                    end
                    modelStructure.start[s] = s0
                    modelStructure.mappedParameters = (;modelStructure.mappedParameters...,  s => s0)
                end
            end


        elseif typeof(v) in [Int64, Float64] || typeof(v) <: Unitful.Quantity || 
                typeof(v) in [Array{Int64,1}, Array{Int64,2}, Array{Float64,1}, Array{Float64,2}] || 
				typeof(v) <: NamedTuple && :class in keys(v) && v.class == :Par ||
				typeof(v) <: NamedTuple && :parameter in keys(v) && v.parameter
            if unitless && !(typeof(v) <: NamedTuple)
                v = ustrip(v)
            end
            modelStructure.parameters[k] = v
            modelStructure.mappedParameters = (;modelStructure.mappedParameters..., k => v)
        elseif (typeof(v) <: NamedTuple && :class in keys(v) && v.class in [:Par, :Var] ||
            typeof(v) <: NamedTuple && :parameter in keys(v) && v.parameter) &&
            :value in keys(v)
            v = v.value
            if typeof(v) in [Expr, Symbol]
                push!(modelStructure.equations, removeBlock(:($k = $v)))
            else
                modelStructure.parameters[k] = v
                modelStructure.mappedParameters = (;modelStructure.mappedParameters..., k => v)
            end
         elseif typeof(v) <: NamedTuple && :class in keys(v) && v.class in [:Var] ||
            typeof(v) <: NamedTuple && :variable in keys(v) && v.variable
            if :init in keys(v)
                x0 = v.init
                if unitless && typeof(x0) != Expr
                    x0 = ustrip(x0)
                end
                modelStructure.init[k] = x0
                modelStructure.mappedParameters = (;modelStructure.mappedParameters..., k => x0)
            end
            if :start in keys(v)
                s0 = v.start
                if unitless && typeof(s0) != Expr
                    s0 = ustrip(s0)
                end
                modelStructure.start[k] = s0
                modelStructure.mappedParameters = (;modelStructure.mappedParameters..., k => s0)
            end
        elseif typeof(v) <: NamedTuple # instantiate
                subModelStructure = ModelStructure()
                flattenModelTuple!(v, subModelStructure, k; unitless, log)
#                println("subModelStructure")
#                printModelStructure(subModelStructure, k)
                mergeModelStructures(modelStructure, subModelStructure, k)
        elseif typeof(v) <:Array && length(v) > 0 && typeof(v[1]) <: NamedTuple # array of instances
            i = 0
            for a in v
                i += 1
                subModelStructure = ModelStructure()
                flattenModelTuple!(a, subModelStructure, k; unitless, log)
                mergeModelStructures(modelStructure, subModelStructure, Symbol(string(k)*"_"*string(i)) )
            end
        elseif isexpr(v, :vect) || isexpr(v, :vcat) || isexpr(v, :hcat)
            arrayEquation = false
            for e in v.args
                if isexpr(e, :(=))
                    if unitless
                        e = removeUnits(e)
                    end
                    push!(modelStructure.equations, removeBlock(e))
                elseif isexpr(e, :tuple)
                    push!(connections, e)
                else
                    arrayEquation = true
                end
            end
            if arrayEquation
                push!(modelStructure.equations, removeBlock(:($k = $(prepend(v, :up)))))
#                 push!(modelStructure.equations, removeBlock(:($k = $v)))
            end                    
        elseif isexpr(v, :(=)) # Single equation
            if unitless
                v = removeUnits(v)
            end
            push!(modelStructure.equations, removeBlock(v))
        elseif v !== nothing # Binding expression
#            println("Binding expression")
#            @show k v typeof(v)

            if unitless
                v = removeUnits(v)
            end
            push!(modelStructure.equations, :($k = $(prepend(v, :up))))
#            push!(modelStructure.equations, :($k = $v))
#            @show modelStructure.equations
        end
    end
#    printModelStructure(modelStructure, "flattened")
    connectEquations = convertConnections!(connections, extendedModel, modelName, log)
    push!(modelStructure.equations, connectEquations...)
end
