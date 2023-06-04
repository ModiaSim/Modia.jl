#=
Recursively instantiate dependent objects and propagate value in hierarchical NamedTuples

* Developer: Hilding Elmqvist, Mogram AB (initial version)
*            Martin Otter, DLR (instantiation of dependent objects added)
* First version: March 2021
* License: MIT (expat)

=#

function getConstructorAsString(path, constructor, parameters):String
    str = ""
    for (key,value) in parameters
        skey = string(key)
        if length(skey) > 0 && skey[1] != '_'
            if str != ""
                str = str*", "
            end
            if typeof(value) == Symbol
                svalue = ":" * string(value)
            elseif typeof(value) == String
                svalue = value
            elseif typeof(value) <: AbstractDict
                svalue = "..."  # Do not show dictionaries, since too much output, especially due to pointers to Object3Ds
            else
                if eltypeOrType(value) <: Number
                    value = ustrip.(value)
                end
                svalue = string(value)
                if length(svalue) > 20
                    svalue = svalue[1:20] * "..."   # Restrict the length of the output
                end
            end
            str = str * skey * "=" * svalue
        end
    end
    if isnothing(path)
        str = string(constructor) * "(" * str * ")"
    else
        str = path * " = " * string(constructor) * "(" * str * ")"
    end
    return str
end


using OrderedCollections: OrderedDict


subst(ex, environment, modelModule) = ex
subst(ex::Symbol, environment, modelModule) = if length(environment) == 0; ex elseif ex in keys(environment[end]); environment[end][ex] else subst(ex, environment[1:end-1], modelModule) end
subst(ex::Expr, environment, modelModule) = Expr(ex.head, [subst(a, environment, modelModule) for a in ex.args]...)
subst(ex::Vector{Symbol}, environment, modelModule) = [subst(a, environment, modelModule) for a in ex]
subst(ex::Vector{Expr}  , environment, modelModule) = [Core.eval(modelModule, subst(a, environment, modelModule)) for a in ex]

#=
function propagate(model, environment=[])
    #println("\n... Propagate: ", model)
    current = OrderedDict()
    for (k,v) in zip(keys(model), model)
        if k in [:_class, :_Type]
            current[k] = v
        elseif typeof(v) <: NamedTuple
            current[k] = propagate(v, vcat(environment, [current]))
        else
            evalv = v
            try
                evalv = eval(subst(v, vcat(environment, [current])))
            catch e
            end

            if typeof(evalv) <: NamedTuple
                current[k] = v
            else
                current[k] = evalv
            end
        end
    end
    return (; current...)
end
=#


appendKey(path, key) = path == "" ? string(key) : path * "." * string(key)


"""
    map = propagateEvaluateAndInstantiate!(partiallyInstantiatedModel::InstantiatedModel; log=false)

Recursively traverse the hierarchical collection `partiallyInstantiatedModel.parameters` and perform the following actions:

- Propagate values.
- Evaluate expressions in the context of `modelModule`.
- Instantiate dependent objects.
- Return the evaluated `parameters` if successfully evaluated, and otherwise
  return nothing, if an error occurred (an error message was printed).
"""
function propagateEvaluateAndInstantiate!(m::InstantiatedModel{FloatType,TimeType}; log=false) where {FloatType,TimeType}
    x_found = fill(false, length(m.equationInfo.x_info))
    map = propagateEvaluateAndInstantiate2!(m, m.parameters, x_found, [], ""; log=log)

    if isnothing(map)
        return nothing
    end

    # Check that all values of x_start are set:
    #x_start_missing = []
    #for (i, found) in enumerate(x_found)
    #    if !found
    #        push!(x_start_missing, eqInfo.x_info[i].x_name)
    #    end
    #end

    # Check that all previous values are set:
    missingInitValues = false
    namesOfMissingValues = ""
    first = true
    for (name,index) in m.previous_dict
        if ismissing(m.previous[index])
            missingInitValues = true
            if first
                first = false
                namesOfMissingValues *= "\n  Variables from previous(..):"
            end
            namesOfMissingValues *= "\n    " * name
        end
    end

    # Check that all pre values are set:
    first = true
    for (name,index) in m.pre_dict
        if ismissing(m.pre[index])
            missingInitValues = true
            if first
                first = false
                namesOfMissingValues *= "\n  Variables from pre(..):"
            end
            namesOfMissingValues *= "\n    " * name
        end
    end

    # Check that all hold values are set:
    first = true
    for (name,index) in m.hold_dict
        if ismissing(m.hold[index])
            missingInitValues = true
            if first
                first = false
                namesOfMissingValues *= "\n  Variables from hold(..):"
            end
            namesOfMissingValues *= "\n    " * name
        end
    end

    if missingInitValues
        printstyled("Model error: ", bold=true, color=:red)
        printstyled("Missing start/init values for variables: ", namesOfMissingValues,
                    bold=true, color=:red)
        print("\n\n")
        return nothing
    end

    if isnothing(map)
        return nothing
    end

    # Resize linear equation systems if dimensions of vector valued tearing variables changed
    resizeLinearEquations!(m, map, m.options.log)

    #if length(x_start_missing) > 0
    #    printstyled("Model error: ", bold=true, color=:red)
    #    printstyled("Missing start/init values for variables: ", x_start_missing,
    #                bold=true, color=:red)
    #    print("\n\n")
    #    return nothing
    #end
    return map
end


"""
    firstName(ex::Expr)

If ex = :(a.b.c.d) -> firstName(ex) = :a
"""
function firstName(ex::Expr)
    if ex.head == :(.)
        if typeof(ex.args[1]) == Expr
            firstName(ex.args[1])
        else
            return ex.args[1]
        end
    end
end

function changeDotToRef(ex)
    if ex.head == :(.)
        ex.head = :ref
        if typeof(ex.args[1]) == Expr
            changeDotToRef(ex.args[1])
        end
    end
    return nothing
end


function propagateEvaluateAndInstantiate2!(m::InstantiatedModel{FloatType,TimeType}, parameters, x_found::Vector{Bool},
                                           environment, path::String; log=false) where {FloatType,TimeType}
    if log
        println("\n 1: !!! instantiate objects of $path: ", parameters)
    end
    current = OrderedDict{Symbol,Any}()   # should be Map()

    # Determine, whether "parameters" has a ":_constructor"  or "_initSegmentFunction" key and handle this specially
    constructor          = nothing
    instantiateFunction  = nothing
    usePath              = false
    useInstantiatedModel = false
    modelModule          = m.modelModule
    eqInfo               = m.equationInfo

    if haskey(parameters, :_constructor)
        # For example: obj = (_class = :Par, _constructor = :(Modia3D.Object3D), _path = true, kwargs...)
        #          or: rev = (_constructor = (_class = :Par, value = :(Modia3D.ModiaRevolute), _path=true), kwargs...)
        v = parameters[:_constructor]
        if typeof(v) <: OrderedDict
            constructor = v[:value]
            if haskey(v, :_path)
                usePath = v[:_path]
            end
            if haskey(v, :_instantiatedmodel)
                useInstantiatedModel = v[:_instantiatedmodel]
            end
        else     
            constructor = v
            if haskey(parameters, :_path)
                usePath = parameters[:_path]
            end
            if haskey(parameters, :_instantiatedModel)
                useInstantiatedModel = parameters[:_instantiatedModel]
            end            
        end
        if log
            if usePath
                println("usePath = true")
            end
            if useInstantiatedModel
                println("useInstantiatedModel = true")
            end            
        end

    elseif haskey(parameters, :_initSegmentFunction)
        # For example: obj = (_initSegmentFunction = Par(functionName = :(instantiateLinearStateSpace!))
        _initSegmentFunction = parameters[:_initSegmentFunction]
        if haskey(_initSegmentFunction, :functionName)
            instantiateFunction = _initSegmentFunction[:functionName]
        else
            @warn "Model $path has key :_initSegmentFunction but its value has no key :functionName"
        end

    elseif haskey(parameters, :value)
        # For example: p1 = (_class = :Var, parameter = true, value = 0.2)
        #          or: p2 = (_class = :Var, parameter = true, value = :(2*p1))
        v = parameters[:value]
        veval = Core.eval(modelModule, subst(v, vcat(environment, [current]), modelModule))
        return veval
    end

    for (k,v) in parameters
        if log
            println(" 2:    ... key = $k, value = $v")
        end
        if k == :_constructor || k == :_buildFunction || k == :_buildOption || k == :_initSegmentFunction || 
           k == :_path || k == :_instantiatedModel || (k == :_class && !isnothing(constructor))
            if log
                println(" 3:    ... key = $k")
            end
            nothing

        elseif !isnothing(constructor) && (k == :value || k == :init || k == :start)
            error("value, init or start keys are not allowed in combination with a _constructor:\n$parameters")

        elseif typeof(v) <: OrderedDict
            if length(v) > 0
                if haskey(v, :_class) && v[:_class] == :Par && haskey(v, :value)
                    # For example: k = (_class = :Par, value = 2.0) -> k = 2.0
                    #          or: k = (_class = :Par, value = :(2*Lx - 3))   -> k = eval( 2*Lx - 3 )
                    #          or: k = (_class = :Par, value = :(bar.frame0)) -> k = ref(bar.frame0)
                    if log
                        println(" 4:    v[:value] = ", v[:value], ", typeof(v[:value]) = ", typeof(v[:value]))
                        #println("        vcat(environment, [current]) = ", vcat(environment, [current]))
                    end
                    subv = subst(v[:value], vcat(environment, [current]), modelModule)
                    if log
                        println(" 5:    _class & value: $k = $subv  # before eval")
                    end
                    if typeof(subv) == Expr && subv.head == :(.)
                        if typeof(firstName(subv)) <: AbstractDict
                            changeDotToRef(subv)
                            if log
                                println(" 5b:    _class & value: $k = $subv  # before eval")
                            end
                        end
                    end
                    current[k] = Core.eval(modelModule, subv)
                    if log
                        println(" 6:                   $k = ", current[k])
                    end
                else
                    if log
                        println(" 7:    ... key = $k, v = $v")
                    end
                    # For example: k = (a = 2.0, b = :(2*Lx))
                    value = propagateEvaluateAndInstantiate2!(m, v, x_found, vcat(environment, [current]), appendKey(path, k); log=log)
                    if log
                        println(" 8:    ... key = $k, value = $value")
                    end
                    if isnothing(value)
                        return nothing
                    end
                    current[k] = value
                end
            end

        else
            if log
                println(" 9:    else: typeof(v) = ", typeof(v))
            end
            subv = subst(v, vcat(environment, [current]), modelModule)
            if log
                println(" 10:          $k = $subv   # before eval")
            end
            if typeof(subv) == Expr && subv.head == :(.)
                if typeof(firstName(subv)) <: AbstractDict
                    changeDotToRef(subv)
                    if log
                        println(" 10b:    _class & value: $k = $subv  # before eval")
                    end
                end
            end

            try
                subv = Core.eval(modelModule, subv)
            catch
                str = getConstructorAsString(path, constructor, parameters)
                printstyled("\nError when instantiating model $(m.modelName):\n$k = $subv\n", bold=true, color=:red)
                Base.rethrow()
            end

            if m.unitless && eltype(subv) <: Number
                # Remove unit
                subv = stripUnit(subv)
            end
            current[k] = subv
            if log
                println(" 11:          $k = ", current[k])
            end

            # Set x_start
            full_key = appendKey(path, k)
            if haskey(eqInfo.x_dict, full_key)
                #if log
                #    println(" 12:              (is stored in x_start)")
                #end
                j = eqInfo.x_dict[full_key]
                xe_info = eqInfo.x_info[j]
                x_value = current[k]
                len = hasParticles(x_value) ? 1 : length(x_value)
                if j <= eqInfo.nx_info_fixedLength && len != xe_info.length
                    printstyled("Model error: ", bold=true, color=:red)
                    printstyled("Length of ", xe_info.x_name, " shall be changed from ",
                                xe_info.length, " to $len\n",
                                "This is not possible because variable has a fixed length.", bold=true, color=:red)
                    return nothing
                end
                x_found[j] = true
                xe_info.startOrInit = deepcopy(x_value)

                # Strip units from x_start
                #if xe_info.length == 1
                #    x_start[xe_info.startIndex] = deepcopy( convert(FloatType, stripUnit(x_value)) )
                #else
                #    ibeg = xe_info.startIndex - 1
                #    for i = 1:xe_info.length
                #        x_start[ibeg+i] = deepcopy( convert(FloatType, stripUnit(x_value[i])) )
                #    end
                #end

            elseif haskey(m.previous_dict, full_key)
                m.previous[ m.previous_dict[full_key] ] = current[k]

            elseif haskey(m.pre_dict, full_key)
                m.pre[ m.pre_dict[full_key] ] = current[k]

            elseif haskey(m.hold_dict, full_key)
                m.hold[ m.hold_dict[full_key] ] = current[k]
            end
        end
    end

    if isnothing(constructor)
        if !isnothing(instantiateFunction)
            # Call: instantiateFunction(model, FloatType, Timetype, buildDict, path)
            # (1) Generate an instance of subModel and store it in buildDict[path]
            # (2) Define subModel states and store them in xxx
            if log
                println(" 13:    +++ Instantiated $path: $instantiateFunction will be called to instantiate sub-model and define varying states\n\n")
            end
            Core.eval(modelModule, :($instantiateFunction($m, $path, $path, $current, log=$log)))
            push!(m.instantiateFunctions, (instantiateFunction, current, path))
        end
        return current
    else
        try
            if usePath
                if useInstantiatedModel
                    obj = Core.eval(modelModule, :(FloatType = $FloatType; $constructor(; path = $path, instantiatedModel = $m, $current...)))                
                else
                    obj = Core.eval(modelModule, :(FloatType = $FloatType; $constructor(; path = $path, $current...)))
                end
            elseif useInstantiatedModel
                obj = Core.eval(modelModule, :(FloatType = $FloatType; $constructor(; instantiatedModel = $m, $current...)))              
            else
                obj = Core.eval(modelModule, :(FloatType = $FloatType; $constructor(; $current...)))
            end
            if log
                println(" 14:    +++ Instantiated $path: typeof(obj) = ", typeof(obj), ", obj = ", obj, "\n\n")
            end
            return obj
        catch
            str = getConstructorAsString(path, constructor, parameters)
            printstyled("\nError in model $(m.modelName) when instantiating\n  $str\n", bold=true, color=:red)
            Base.rethrow()
        end
    end
end


"""
    splittedPath = spitPlath(path::Union{Symbol, Expr, Nothing})::Vector{Symbol}

# Examples
```
splitPath(nothing)     # = Symbol[]
splitPath(:a)          # = Symbol[:a]
splitPath(:(a.b.c.d))  # = Symbol[:a, :b, :c, :d]
```
"""
function splitPath(path::Union{Symbol, Expr, Nothing})::Vector{Symbol}
    splittedPath = Symbol[]
    if typeof(path) == Symbol
        push!(splittedPath, path)
    elseif typeof(path) == Expr
        while typeof(path.args[1]) == Expr
            pushfirst!(splittedPath, path.args[2].value)
            path = path.args[1]
        end
        pushfirst!(splittedPath, path.args[2].value)
        pushfirst!(splittedPath, path.args[1])
    end
    return splittedPath
end


"""
    modelOfPath = getModelFromSplittedPath(model, splittedPath::Vector{Symbol})

Return reference to the sub-model characterized by `splittedPath`.

# Examples

```
mymodel = Model(a = Model(b = Model(c = Model)))
model_b = getModelFromSplittedPath(mymodel, Symbol["a", "b"])  # = mymodel[:a][:b]
```
"""
function getModelFromSplittedPath(model, splittedPath::Vector{Symbol})
    for name in splittedPath
        model = model[name]
    end
    return model
end