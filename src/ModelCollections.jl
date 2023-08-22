#=
Handles models and variables defined as dictionaries.

* Developer: Hilding Elmqvist, Mogram AB
* First version: June 2021
* License: MIT (expat)

=#

export mergeModels, recursiveMerge, redeclare, outer, showModel, @showModel, drawModel, Model, Map, Par, Var, setLogMerge,
    constant, parameter, input, output, potential, flow, interval, @info_str, Boolean, Integ, @define, stringifyDefinition, Lookup,
    resetLibrary, printLibrary, makePlacement

using Base.Meta: isexpr
using OrderedCollections: OrderedDict
using Unitful
using ModiaBase.Symbolic: removeBlock, prepend
#using Profile

global logMerge = false

function setLogMerge(val)
    global logMerge
    logMerge = val
end

Lookup(model) = (:Lookup, model)

function stringifyDefinition(v)
    if typeof(v) in [Symbol, Expr]
        v = removeBlock(v)
        v = ":(" * string(v) * ")"
    end
    return v
end

"Quote an expression; don't quote if the meaning in the AST is the same anyway."
selquote(x) = x
selquote(sym::Symbol) = Meta.quot(sym)
function selquote(ex::Expr)
    if isexpr(ex, :kw) || isexpr(ex, :call)
        Expr(ex.head, [ex.args[1], [selquote(e) for e in ex.args[2:end]]...]...) 
    else
        Expr(ex.head, [selquote(e) for e in ex.args]...)
    end
end

macro define(modelDef)
    modelName = modelDef.args[1]
    println(modelName)
    model = modelDef.args[2]
    model = selquote(model)
    dump(model)
    res = :($modelName = Meta.quot($model))
    return esc(res)
end

#Base.:<(m::AbstractDict, n::AbstractDict) = m # OrderedDict()

# Convert hierchical modifiers and dot notation to merge operator
mergify(ex, model, kwvalue=false, level=1) = ex
function mergify(ex::Expr, model, kwvalue=true, level=1)
    println("\nmergify: ")
    level += 1
    @show ex model
#    showModel(ex)
    res = nothing
#    dump(ex)
    if ex.args[1] == :equations || ex.args[1] == :*  || ex.head == :ncat # TODO: refine
        res = ex
    elseif isexpr(ex, :call) && !(ex.args[1] in [:Model, :Map, :Par, :Var, :|, :connect])
#        println("Instantiation")
        mapArguments = [mergify(e, model, ex.head == :kw, level) for e in ex.args[2:end]]
        if length(mapArguments) > 0
            maps = foldl( (x,y) -> :($x | Map($y)), mapArguments[2:end], init=:(Map($(mapArguments[1]))))
            println("\nconverted: ", showModel(mapArguments))
            println("  to: ", showModel(maps))
        end
        if kwvalue 
            redec = level > 4  # TODO: Better criteria
            if redec
                if length(mapArguments) > 0
                    res = :(redeclare | $(ex.args[1]) | $maps) # Map($(mapArguments...)) )
                else
                    res = :(redeclare | $(ex.args[1]) )
                end
            else
                if length(mapArguments) > 0
#                    @show mapArguments
                    res = :($(ex.args[1]) | $maps) # Map($(mapArguments...)) )
                else
#                    @show ex.args[1] 
#                    dump(model)
                    modelDict = OrderedDict([kv.args[1] => kv.args[2] for kv in model.args[2:end]]...)
#                    @show modelDict
                    compType = ex.args[1]
#                    @show compType
                    if compType in keys(modelDict)
#=
                        newType = modelDict[compType]
                        @show newType
                        dump(newType)
                        compType = :($(newType.args[2]))
=#
#                        res =  :((:Lookup, $(Meta.quot(compType))))
                        res =  :(Lookup($(Meta.quot(compType))))
#                        println("mergified: ") 
#                        showModel(res)                    
                        return res
                    end
                    res = :($(compType) )
                end
            end
        else
            func = ex.args[1]
#            @show ex func
            if func != :(=>)
                res = Expr(:kw, func, :(Map($(mapArguments...))) )  
            else
                dump(ex)
#                :($(ex.args[2].args[1]) = 0)  
                res = Expr(:kw, ex.args[2].args[1], Expr(:kw, ex.args[2].args[1], 0) )  
            end
        end
    elseif ex.head == :kw && typeof(ex.args[1]) == Expr && length(ex.args) == 2 && ex.args[1].head == :.
#        println("Dot notation")
        sequence = []
        unpackPath(ex.args[1], sequence)
        nest = foldr( (x,y) -> :(Map($x=$y)), sequence[2:end], init=ex.args[2])
        res = Expr(:kw, sequence[1], nest)
#        println("converted dot notation: ", ex)
#        println("  to: ", res)
#        Expr(:call, :(=>), [mergify(e, model, ex.head == :kw, level) for e in ex.args]...)
    else
#        println("Model, Map, Par, Var, |, connect or :kw: ", ex.args[1])
        res = Expr(ex.head, [mergify(e, model, ex.head == :kw, level) for e in ex.args]...)
    end
#    println("mergified: ") 
#    showModel(res)
    res
end

newLibrary = OrderedDict{Symbol, Any}()

function makeVar(type, type_prefix = "")
    v = OrderedDict{Symbol, Any}()
    v[:type_specifier] = type
    if type_prefix != ""
        v[:type_prefix] = type_prefix
    end
    v
end


MakeMap(; kwargs...) = (; kwargs...)

function makeComp(constructor) 
    @show constructor
    dump(constructor)
    class_prefixes = []
    if constructor == :potential
        makeVar("Real")
    elseif constructor == :flow
        makeVar("Real", ["flow"])
    elseif typeof(constructor) == Float64
        v = makeVar("Real")
        v[:type_prefix] = ["parameter"]
        v[:binding] = constructor
        v
    elseif typeof(constructor) == String
        v = makeVar("String")
        v[:type_prefix] = ["parameter"]
        v[:binding] = constructor
        v
    elseif typeof(constructor) == Expr && constructor.args[1] == :* && constructor.args[3].head == :macrocall
        v = makeVar("Real")
        v[:type_prefix] = ["parameter"]
        v[:binding] = constructor.args[2]
        v[:class_modification] = OrderedDict(:unit => OrderedDict(:binding => constructor.args[3].args[3]))
        v
    elseif hasproperty(constructor, :head) && constructor.head == :call && constructor.args[1] == :|
        makeVar(constructor.args[3])
    elseif hasproperty(constructor, :head) && constructor.head == :call && constructor.args[1] == :Var
        # class_modification!!!
        OrderedDict(constructor.args[2].args[1] => constructor.args[2].args[2])
    elseif hasproperty(constructor, :head) && constructor.head == :call && constructor.args[2].head == :kw && constructor.args[2].args[2].value == :left
        v = makeVar(constructor.args[1])
        v[:annotation] = MakeMap(
            annotation = MakeMap(
                Placement = MakeMap(
                    class_modification = MakeMap(
                        transformation = MakeMap(
                            class_modification = MakeMap(
                                extent = MakeMap(
                                binding = "[[-110,-10],[-90,10]]"
                                )
                            )
                        )
                    )
                    )
                )
            )
        v
    elseif hasproperty(constructor, :head) && constructor.head == :call && constructor.args[2].head == :kw && constructor.args[2].args[2].value == :right
        v = makeVar(constructor.args[1])
        v[:annotation] = MakeMap(
            annotation = MakeMap(
                Placement = MakeMap(
                    class_modification = MakeMap(
                        transformation = MakeMap(
                            class_modification = MakeMap(
                                extent = MakeMap(
                                binding = "[[110,-10],[90,10]]"
                                )
                            )
                        )
                    )
                    )
                )
            )
        v
    else
        makeVar(constructor)
    end
end

function makeAST(JuliaAST)
    @show JuliaAST
    ast = OrderedDict{Symbol, Any}()
    class_prefixes = []
    elements = OrderedDict{Symbol, Any}()
    if hasproperty(JuliaAST, :head) && JuliaAST.head == :call && JuliaAST.args[1] == :|
        println("\n\n\nHEJSAN")
        @show JuliaAST.args[2]
        ast = makeAST(JuliaAST.args[3])
        ast[:extends] = OrderedDict(:type_specifier => JuliaAST.args[2])
    else
        for a in JuliaAST.args[2:end]
            @show a
            n = a.args[1]
            if n != :equations
                constructor = a.args[2]
                elements[n] = makeComp(constructor) 
                if :type_prefix in keys(elements[n]) && elements[n][:type_prefix] == ["flow"]
                    class_prefixes = ["connector"]
                end
            end
            @show elements
        end
        if class_prefixes != []
            ast[:class_prefixes] = class_prefixes
        end
        ast[:elements] = elements
    end
    ast
end

function storeModelAST(modelDef)
    dump(modelDef)
    modelName = modelDef.args[1]
    JuliaAST = modelDef.args[2]
    newLibrary[modelName] = makeAST(JuliaAST)
end

function resetLibrary()
    global newLibrary
    newLibrary = OrderedDict{Symbol, Any}()
end

function printLibrary()
    println("newLibrary")
    println(JSON.json(newLibrary, 2))
end

macro defineNEW(modelDef)
#    storeModelAST(modelDef)
    println()
    modelName = modelDef.args[1]
    println("Transforming: ", modelName)
    model = modelDef.args[2]
    println("Original model:")
    print("@define ", modelName, " = ")
    showModel(model)
    println()
#    dump(model, maxdepth=14)
#    model = selquote(model)
    model = mergify(model, model)
    println("Modified model:")
    print(modelName, " = ")
    showModel(model)
    println()
#    dump(model, maxdepth=14)
    res = :($modelName = $model) # Meta.quot($model))
    return esc(res)
end

function substituteGenerics(model)
    for (k,v) in model
        if isCollection(v)
            model[k] = substituteGenerics(v)
        elseif typeof(v) <: Tuple && v[1] == :Lookup
            # (:_lookup, name [, modifier])
            name = v[2]
            # lookup name and set componenttype
#            @show model[name]
            componentType = deepcopy(model[name])
#            @show name k componentType
#            println("Redeclare $k = ")
#            showModel(componentType)
            componentType[:_class] = :Model
            delete!(componentType, :value)
#            model[name] = OrderedDict()
            model[name][:_class] = :Par
            if length(v) == 3
                model[k] = componentType | v[3]
            else
                model[k] = componentType
            end
#            @show name model[k]
        end
    end
    model
end

# Remove generic parameters after all lookups have been performed
function removeGenerics(model)
    for (k,v) in model
        if isCollection(v)
            if v[:_class] == :Par
                model[k] = OrderedDict() # Replace by empty OrderedCollection
#                delete!(model, k)
            else
                model[k] = removeGenerics(v)
            end
        end
    end
    model
end

function showModel(m::Expr)
    if isexpr(m, :call) && m.args[1] == :Model
        println("Model(")
        for kw in m.args[2:end]
            println("  ", kw.args[1], " = ", kw.args[2])
        end
        println(")")
    elseif isexpr(m, :kw)
        println("  ", m.args[1], " = ", m.args[2])
    else
        println(m)
    end
end

isCollection(v) = (typeof(v) <: AbstractDict) && :_class in keys(v)

showModel(m, level=0) = println(m)
function showModel(m::OrderedDict, level=0)
#    @show m typeof(m)
    level += 1
    if isCollection(m)
        print(m[:_class])
    else
        print("OrderedDict")
    end
    println("(")
    for (k, v) in m
        if isCollection(v) 
            print("  "^level, k, " = ")
            showModel(v, level)
        elseif k == :equations
            println("  "^level, "equations = :[")  
            for e in v.args
                println("  "^(level+1), e)
            end
            println("  "^level, "]")         
        elseif typeof(v) <: AbstractArray
            if eltype(v) <: Number
                println("  "^level, k, " = ", v)
            else
                println("  "^level, k, " = [")  
                for e in v
                    print("  "^(level+1))
                    showModel(e, level+1)
                end
                println("  "^level, "]")
            end
        elseif k != :_class
            #println("  "^level, k, " = ", stringifyDefinition(v), ",")
            print("  "^level, k, " = ")
            showModel(v,level)
        end
    end
    println("  "^(level-1), "),")
end

macro showModel(model)
    modelName = string(model)
    mod = :( print($modelName, " = "); showModel($model); println() )
    return esc(mod)
end

id = 1

function drawModel(name, m, level=0)
    this = Dict()
    children = []
    ports = []
    edges = []
    level += 1
    global id
    if typeof(m) <: AbstractDict && haskey(m, :_class)
#        print(m._class)
    end
    for (k, v) in zip(keys(m), m)
        if typeof(v) <: AbstractDict && level <= 1
            child = drawModel(k, v, level)
            push!(children, child)
        elseif k in [:_class, :equations, :input, :output, :potential, :flow, :info]
        elseif k == :connect
            for conn in v.args
                conn = conn.args
                source = conn[1]
                target = conn[2:end]
                id += 1
                edge = (; id, sources = [string(source)], targets=[string(target[1])])
                push!(edges, edge)
            end
        elseif level <= 2
            port = (; id=string(name) * "." * string(k), height=10, width=10)
            push!(ports, port)
        end
    end

    if length(children) > 0 && length(ports) > 0
        this = (;id=name, height=100, width=100, labels=[(;text=name)], ports, children, edges)
    elseif length(children) > 0 
        this = (;id=name, height=100, width=100, labels=[(;text=name)], children, edges)
    elseif length(ports) > 0 
        this = (;id=name, height=100, width=100, labels=[(;text=name)], ports)
    else
        this = (;id=name, height=100, width=100, labels=[(;text=name)])
    end
    return this
end


# ---------------------------------------------------------

function mergeModels(m1::AbstractDict, m2::AbstractDict, env=Symbol())
#    println("mergeModels")
#    @show m1 m2
    result = deepcopy(m1) 
    for (k,v) in m2 
        if k == :_class
        elseif typeof(v) <: AbstractDict
            if k in keys(result) && ! (:_redeclare in keys(v))
                if logMerge; print("In $k: ") end
#                @show result[k] v k
                if typeof(result[k]) <: AbstractDict
                    m = mergeModels(result[k], v, k)
                    result[k] = m
                else
#                    result[k] = v
                    result[k] = result[k] | v # Fails with @define on LinearSystem
                end
            elseif :_redeclare in keys(v)
                if logMerge; println("Redeclaring: $k = $v") end
                result[k] = v
            elseif nothing in values(v) # TODO: Refine

            else
                if !(:_redeclare in keys(result))
                    if logMerge; print("Adding: $k = "); showModel(v, 2) end
                end
                result[k] = v
            end
        elseif v === nothing
            if logMerge; println("Deleting: $k") end
            delete!(result, k)
        elseif k in keys(result) && k == :equations
            equa = copy(result[k])
            push!(equa.args, v.args...)
            result[k] = equa
            if logMerge
                println("Adding equations: ", v)
            end
        else
            if logMerge
                if k in keys(result)
                    if result[k] != v
                        println("Changing: $k = $(stringifyDefinition(result[k])) to $k = $(stringifyDefinition(v))")
                    end
                elseif !(:_redeclare in keys(result))
                    println("Adding: $k = $(stringifyDefinition(v))")
                end
            end
            result[k] = v
        end
    end
#    delete!(result, :_class)
    return result 
end

function newCollection(kwargs, kind)
    m = OrderedDict{Symbol, Any}(kwargs)
    m[:_class] = kind
    m 
end

Model(; kwargs...) = newCollection(kwargs, :Model)

Map(; kwargs...) = newCollection(kwargs, :Map)

Par(; kwargs...) = newCollection(kwargs, :Par)

Var(;kwargs...) = newCollection(kwargs, :Var)

par = Par()
var = Var()
type = Map()

function Var(values...; kwargs...)
    res = newCollection(kwargs, :Var) 
    for v in values
        res = res | v
    end
    res
end

function Par(value; kwargs...)
    res = newCollection(kwargs, :Par)
    res[:value] = value
    res
end

Integ(;kwargs...) = Var(kwargs)

Boolean(;kwargs...) = Var(kwargs)

redeclare = Model( _redeclare = true)
redefine = Par( _redeclare = true)
outer = Var( _outer = true)

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

Base.:|(m::AbstractDict, n::AbstractDict) =  mergeModels(m, n) 
Base.:|(m::Symbol, n::AbstractDict) =  mergeModels(eval(m), n) 
Base.:|(m, n) = if !(typeof(n) <: AbstractDict); mergeModels(m, Var(value=n)) else mergeModels(n, Var(value=m)) end
Base.:|(m::Tuple, n::AbstractDict) =  (m[1],  m[2], n)

#mergeModels(x, ::Nothing) = x
#mergeModels(x, y) = y
#recursiveMerge(x::Expr, y::Expr) = begin dump(x); dump(y); Expr(x.head, x.args..., y.args...) end
#mergeModels(x::Expr, y::Tuple) = begin x = copy(x); xargs = x.args; xargs[y[2]] = y[3]; Expr(x.head, xargs...) end


# -------------------------------------------------------------------------------

# Handling of connections

"""
Unpacks a path to a sequence.
Example: :(a.b.c) => [:a, :b, :c]
"""
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

isVar(v, kind) = isCollection(v) && v[:_class] == :Var && kind in keys(v) && v[kind]

"""
Collects potentials, flows, inputs and outputs from a Model.
"""
function collectConnector(name, model)
    potentials = []
    flows = []
    inputs = OrderedDict{Any,Any}()
    outputs = OrderedDict{Any,Any}()
    potentialPotentials = []
    for (k,v) in model
        if isVar(v, :potential)
            push!(potentials, k)
        elseif isVar(v, :flow)
            push!(flows, k)
            push!(potentials, potentialPotentials...)
        elseif isVar(v, :input)
            inputs[prepend(k, name)] = v 
        elseif isVar(v, :output)
            outputs[prepend(k, name)] = v
        elseif isCollection(v) && v[:_class] == :Var
#            println("Non potential or flow variable: $k = $v")
            push!(potentialPotentials, k)
        elseif isCollection(v) && v[:_class] == :Model
            subPotentials, subFlows = collectConnector(k, v)
            push!(potentials, [prepend(p, k) for p in subPotentials]...)
            push!(flows, [prepend(f, k) for f in subFlows]...)
        end
    end
    return potentials, flows, inputs, outputs
end

"""
Merge connections referencing the same connector into the same tuple (similar as Modelica semantics of connections).
For examples, see below.
"""
function mergeConnections!(connections)
    for i in 1:length(connections)
        for j in 1:i-1
            con1 = connections[i]
            con2 = connections[j]
            if length(con1.args) > 0 && length(con2.args) > 0 && length(intersect(Set(con1.args), Set(con2.args))) > 0
                connections[i] = Expr(:tuple, union(Set(con1.args), Set(con2.args))...)
                connections[j] = Expr(:tuple) # Empty tuple
            end
        end
    end
end

function testMergeConnections()
    connections = [
        :((a, b))
        :((b, c))
        :((m1.b, m2.a))
        :((m1.b, m3.a))
        :((r,s,t))
        :((s,t,u))
    ]
    mergeConnections!(connections)
    @assert connections == Expr[:(()), :((a, b, c)), :(()), :((m1.b, m3.a, m2.a)), :(()), :((s, u, t, r))]
end

"""
Convert connections to equations
"""
function convertConnections!(connections, model, modelName, logging=false)
#    println("\nconvertConnections")
#    showModel(model)
    mergeConnections!(connections)

    connectEquations = []
    alreadyConnected = []
    for i in 1:length(connections)
        c = connections[i]
        if c.head == :tuple
            connected = c.args

            inflows = []
            outflows = []
            first = true
            signalFlow1 = nothing
            connectedOutput = nothing
            potentials1 = nothing
            fullPotentials1 = []
            flows1 = []
            for con in connected
                if con in alreadyConnected
                    error("Already connected: $con, found in connection: $connected")
                end
                push!(alreadyConnected, con)
#                @show connected con

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
                    if :input in keys(mod) && mod[:input] || :output in keys(mod) && mod[:output] || :_class in keys(mod) && mod[:_class] == :Var
                        signalFlow = con
                        if signalFlow1 !== nothing
                            push!(connectEquations, :($signalFlow1 = $signalFlow))
                        end
                        signalFlow1 = signalFlow

                        if :output in keys(mod) && mod[:output]
                            if connectedOutput != nothing # TODO: refine logic concerning inner and outer inputs and outputs
#                                error("It is not allowed to connect two outputs: ", connectedOutput, ", ", con)
                            end
                            connectedOutput = con
                        end
                    elseif :_class in keys(mod) && mod[:_class] == :Var
                        if length(fullPotentials1) > 0
                            push!(connectEquations, :($(fullPotentials1[1]) = $con))
                        end
                        push!(fullPotentials1, con)
                    else
                        potentials, flows = collectConnector("", mod)

                        # Deprecated
                        if :potentials in keys(mod)
                            potentials = vcat(potentials, mod.potentials.args)
                        end
                        if :flows in keys(mod)
                            flows = vcat(flows, mod.flows.args)
                        end
                        
                        if potentials1 != nothing && potentials != potentials1
                            error("Not compatible potential variables: $potentials1 != $potentials, found in connection: $connected")
                        end
                        fullPotentials = []
                        for i in 1:length(potentials)
                            p = potentials[i]
                            potential = prepend(p, con)
                            push!(fullPotentials, potential)
                            if potentials1 != nothing
                                push!(connectEquations, :($(fullPotentials1[i]) = $potential))
                            end
                        end
                        potentials1 = potentials 
                        fullPotentials1 = fullPotentials

                        if first 
                            for i in 1:length(flows)
                                push!(inflows, 0)
                                push!(outflows, 0)
                            end
                            first = false
                        end

                        if length(flows1) > 0 && flows != flows1
                            error("Not compatible flow variables: $flows1 != $flows, found in connection: $connected")
                        end
                        for i in 1:length(flows)
                            f = flows[i]
                            flowVar = prepend(f, con)
                            if length(sequence) == 1
                                if inflows[i] == 0
                                    inflows[i] = flowVar
                                else
                                    inflows[i] = :($(inflows[i]) + $flowVar)
                                end
                            else
                                if outflows[i] == 0
                                    outflows[i] = flowVar
                                else
                                    outflows[i] = :($(outflows[i]) + $flowVar)
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
            for i in 1:length(inflows)
                if inflows[i] != 0 || outflows[i] != 0
                    push!(connectEquations, :($(inflows[i]) = $(outflows[i])))
                end
            end
        end
    end

    if length(connectEquations) > 0 && logging
#        testMergeConnections()

        println("Connect equations in $modelName:")
        for e in connectEquations
            println("  ", e)
        end
    end
    return connectEquations
end

# ----------------------------------------------------------------------------

# Flattening

function flattenModelTuple!(model, modelStructure, modelName, to; unitless = false, log=false)
#    @show model
    connections = []
    extendedModel = deepcopy(model)

    for (k,v) in model 
        if typeof(modelName) == String
            subMod = k
        else
            subMod = :($modelName.$k)
        end
        if k in [:_class, :_redeclare]
            # skip them
        elseif typeof(v) <: Number || typeof(v) <: Unitful.Quantity || 
                typeof(v) in [Array{Int64,1}, Array{Int64,2}, Array{Float64,1}, Array{Float64,2}] || 
				isCollection(v) && v[:_class] == :Par 
#				typeof(v) <: AbstractDict && :parameter in keys(v) && v.parameter
            # p = 5.0
            if unitless && !isCollection(v)
                v = ustrip.(v)
            end
            modelStructure.parameters[subMod] = v
            modelStructure.mappedParameters[k] = v
        elseif (isCollection(v) && v[:_class] in [:Par, :Var] ||
                isCollection(v) && :parameter in keys(v) && v[:parameter]) &&
                :value in keys(v)
            v = v[:value]
            if typeof(v) in [Expr, Symbol]
                push!(modelStructure.equations, removeBlock(:($k = $v)))
            else
                modelStructure.parameters[subMod] = v
                modelStructure.mappedParameters[k] = v
            end
         elseif isCollection(v) && v[:_class] in [:Var] ||
            typeof(v) <: AbstractDict && :variable in keys(v) && v[:variable]
            if :init in keys(v) && !isnothing(v[:init])
                x0 = v[:init]
                if unitless && typeof(x0) != Expr
                    x0 = ustrip.(x0)
                end
                modelStructure.init[subMod] = x0
                modelStructure.mappedParameters[k] = x0
            end
            if :start in keys(v) && !isnothing(v[:start])
                s0 = v[:start]
                if unitless && typeof(s0) != Expr
                    s0 = ustrip.(s0)
                end
                modelStructure.start[subMod] = s0
                modelStructure.mappedParameters[k] = s0 
            end
            if :input in keys(v) && v[:input]
                modelStructure.inputs[subMod] = v
            end
            if :output in keys(v) && v[:output]
                modelStructure.outputs[subMod] = v
            end
            if :_outer in keys(v) && v[:_outer]
                push!(modelStructure.equations, :($k = $(prepend(k, :up))))
            end
            if :hideResult in keys(v) && v[:hideResult] == true
                push!(modelStructure.hideResults, subMod)
            end
        elseif isCollection(v) # || typeof(v) == Symbol # instantiate
                if typeof(v) == Symbol
                    v = eval(eval(v))
                end
                subModelStructure = ModelStructure()
                m = :($modelName.$k)
                if typeof(modelName) == String
                    m = k
                end
                flattenModelTuple!(v, subModelStructure, m, to; unitless, log)
#=
                println("subModelStructure")
                @show subModelStructure
                printModelStructure(subModelStructure, label=k)
=#
            @timeit to "merge model structures" mergeModelStructures(modelStructure, subModelStructure, k)
        elseif typeof(v) <:Array && length(v) > 0 && (typeof(v[1]) <: AbstractDict || typeof(v[1]) <: AbstractDict) # array of instances
            i = 0
            for a in v
                i += 1
                subModelStructure = ModelStructure()
                name = Symbol(string(k)*"_"*string(i))
                flattenModelTuple!(a, subModelStructure, name, to; unitless, log)
                @timeit to "merge model structures" mergeModelStructures(modelStructure, subModelStructure, name )
            end
        elseif isexpr(v, :vect) || isexpr(v, :vcat) || isexpr(v, :hcat)
            arrayEquation = false
            for e in v.args
                if k == :equations && isexpr(e, :(=)) # Allow expressions without result
                    if unitless
                        e = removeUnits(e)
                    end
                    push!(modelStructure.equations, removeBlock(e))
                elseif isexpr(e, :tuple)
                    push!(connections, e)
                elseif isexpr(e, :call) && e.args[1] == :connect
                    con =  Expr(:tuple, e.args[2:end]... )
                    push!(connections, con)
                else
                    arrayEquation = true
                end
            end
            if arrayEquation
                push!(modelStructure.equations, removeBlock(:($k = $(prepend(v, :up)))))  # Needed for TestStateSpace.SecondOrder1
#                push!(modelStructure.equations, removeBlock(:($k = $v)))   # To handle assert
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
#            push!(modelStructure.equations, :($k = $v))  # To fix Modelica_Electrical_Analog_Interfaces_ConditionalHeatPort
#            @show modelStructure.equations
        end
    end

    @timeit to "convert connections" connectEquations = convertConnections!(connections, extendedModel, modelName, log)
    push!(modelStructure.equations, connectEquations...)

    printModelStructure(modelStructure, "flattened", log=false)
end
