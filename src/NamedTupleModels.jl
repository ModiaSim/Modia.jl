#=
Handles models defined as named tuples.

* Developer: Hilding Elmqvist, Mogram AB
* First version: January 2021
* License: MIT (expat)

=#
export mergeModels, recursiveMerge, Redeclare, showModel, @showModel, Model, Map, setLogMerge

using Base.Meta: isexpr
using DataStructures: OrderedDict
using Unitful
using ModiaBase.Symbolic: removeBlock, prepend


function showModel(m, level=0)
    println("(")
    level += 1
    for (k, v) in zip(keys(m), m)
        if typeof(v) <: NamedTuple
            print("  "^level, k, " = ")
            showModel(v, level)
        else
            println("  "^level, k, " = ", removeBlock(v), ",")
        end
    end
    println("  "^(level-1), "),")
end

macro showModel(model)
    modelName = string(model)
    mod = :( print($modelName, " = Model"); showModel($model); println() )
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
            elseif nothing in values(v)  # TODO: Refine

            else
                if !(:_redeclare in keys(mergedModels))
                    if logMerge; println("Adding: $k = $v") end
                end
                mergedModels[k] = v
            end
        elseif v === nothing
            if logMerge; println("Deleting: $k") end
            delete!(mergedModels, k)
        else
            if logMerge
                if k in keys(mergedModels)
                    println("Changing: $k = $(mergedModels[k]) to $k = $v")
                elseif !(:_redeclare in keys(mergedModels))
                    println("Adding: $k = $v")
                end
            end
            mergedModels[k] = v
        end
    end
    return (; mergedModels...) # Transform OrderedDict to named tuple
end

Model(; kwargs...) = (; kwargs...)

Map(; kwargs...) = (; kwargs...) # OrderedDict{Symbol,Any}(kwargs)

Base.:∪(m::NamedTuple, n::NamedTuple) =  mergeModels(m, n)
Base.:|(m::NamedTuple, n::NamedTuple) =  mergeModels(m, n)

Redeclare = ( _redeclare = true, )

recursiveMerge(x, ::Nothing) = x
recursiveMerge(x, y) = y
recursiveMerge(x::Expr, y::Expr) = begin dump(x); dump(y); Expr(x.head, x.args..., y.args...) end
recursiveMerge(x::Expr, y::Tuple) = begin x = copy(x); xargs = x.args; xargs[y[2]] = y[3]; Expr(x.head, xargs...) end

function recursiveMerge(nt1::NamedTuple, nt2::NamedTuple)
    all_keys = union(keys(nt1), keys(nt2))
    gen = Base.Generator(all_keys) do key
        v1 = get(nt1, key, nothing)
        v2 = get(nt2, key, nothing)
        key => recursiveMerge(v1, v2)
    end
    return Map(; gen...)
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
                    potentials = mod.potentials.args
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

                    flows = mod.flows.args
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
                else
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
        if k in [:inputs, :outputs, :potentials, :flows]

        elseif k == :init
            for (x,x0) in zip(keys(v), v)
                if unitless && typeof(x0) != Expr
                    x0 = ustrip(x0)
                end
                modelStructure.init[x] = x0
                modelStructure.mappedParameters = (;modelStructure.mappedParameters..., x => x0)
            end
        elseif k == :start
            for (s,s0) in zip(keys(v), v)
                if unitless
                    s0 = ustrip(s0)
                end
                modelStructure.start[s] = s0
                modelStructure.mappedParameters = (;modelStructure.mappedParameters...,  s => s0)
            end
        elseif typeof(v) in [Int64, Float64] || typeof(v) <: Unitful.Quantity || typeof(v) in [Array{Float64,1}, Array{Float64,2}]
            if unitless
                v = ustrip(v)
            end
            modelStructure.parameters[k] = v
            modelStructure.mappedParameters = (;modelStructure.mappedParameters..., k => v)
        elseif typeof(v) <: NamedTuple # instantiate
                subModelStructure = ModelStructure()
                flattenModelTuple!(v, subModelStructure, k; unitless, log)
#                println("subModelStructure")
#                printModelStructure(subModelStructure, k)
                mergeModelStructures(modelStructure, subModelStructure, k)
        elseif typeof(v) <:Array && length(v) > 0 && typeof(v[1]) <: NamedTuple
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
            end                    
        elseif isexpr(v, :(=)) # Single equation
            if unitless
                v = removeUnits(v)
            end
            push!(modelStructure.equations, removeBlock(v))
        elseif v !== nothing # Binding equation
            if unitless
                v = removeUnits(v)
            end
            push!(modelStructure.equations, :($k = $(prepend(v, :up))))
#            @show modelStructure.equations
        end
    end
#    printModelStructure(modelStructure, "flattened")
#    @show extendedModel
    connectEquations = convertConnections!(connections, extendedModel, modelName, log)
    push!(modelStructure.equations, connectEquations...)
end
