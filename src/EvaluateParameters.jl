#=
Recursively instantiate dependent objects and propagate value in hierarchical NamedTuples

* Developer: Hilding Elmqvist, Mogram AB (initial version)
*            Martin Otter, DLR (instantiation of dependent objects added)
* First version: March 2021
* License: MIT (expat)

=#


using DataStructures: OrderedDict


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
        if k in [:class, :_Type]
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
    propagateEvaluateAndInstantiate(modelModule::Module, model::NamedTuple)
    
Recursively traverse the hierarchical NamedTuple `model` and perform the following actions:

- Propagate values.
- Evaluate expressions in the context of `modelModule`.
- Instantiate dependent objects.
- Return the evaluated `model` as NamedTuple.
"""
function propagateEvaluateAndInstantiate(modelModule, model, environment=[], path=""; log=false)
    if log
        println("\n!!! instantiate objects of $path: ", model)
    end
    current = OrderedDict()
    
    # Determine, whether the "model" has a ":_constructor" key and handle this specially
    constructor = nothing
    usePath     = false
    if haskey(model, :_constructor)
        # For example: obj = (class = :Par, _constructor = :(Modia3D.Object3D), _path = true, kwargs...)
        #          or: rev = (_constructor = (class = :Par, value = :(Modia3D.ModiaRevolute), _path=true), kwargs...)    
        v = model[:_constructor]
        if typeof(v) <: NamedTuple
            constructor = v[:value]
            if haskey(v, :_path)
                usePath = v[:_path]
            end
        else
            constructor = v
            if haskey(model, :_path)
                usePath = model[:_path]
            end
        end
        
    elseif haskey(model, :value)
        # For example: p1 = (class = :Var, parameter = true, value = 0.2)
        #          or: p2 = (class = :Var, parameter = true, value = :(2*p1))
        v = model[:value]
        veval = Core.eval(modelModule, subst(v, vcat(environment, [current]), modelModule))
        return veval
    end
    
    for (k,v) in zip(keys(model), model)
        if log
            println("    ... key = $k, value = $v")
        end
        if k == :_constructor || k == :_path || (k == :class && !isnothing(constructor))
            nothing
            
        elseif !isnothing(constructor) && (k == :value || k == :init || k == :start)
            error("value, init or start keys are not allowed in combination with a _constructor:\n$model")
            
        elseif typeof(v) <: NamedTuple    
            if haskey(v, :class) && v[:class] == :Par && haskey(v, :value)
                # For example: k = (class = :Par, value = 2.0) -> k = 2.0
                #          or: k = (class = :Par, value = :(2*Lx - 3))   -> k = eval( 2*Lx - 3 )   
                #          or: k = (class = :Par, value = :(bar.frame0)) -> k = ref(bar.frame0)
                subv = subst(v[:value], vcat(environment, [current]), modelModule)
                if log
                    println("    class & value: $k = $subv  # before eval")
                end
                current[k] = Core.eval(modelModule, subv)
                if log
                    println("                   $k = ", current[k])
                end
            else
                # For example: k = (a = 2.0, b = :(2*Lx))
                current[k] = propagateEvaluateAndInstantiate(modelModule, v, vcat(environment, [current]), appendKey(path, k); log=log)     
            end
            
        else
            if log
                println("    else: typeof(v) = ", typeof(v))
            end
            subv = subst(v, vcat(environment, [current]), modelModule)
            if log
                println("          $k = $subv   # before eval")
            end
            current[k] = Core.eval(modelModule, subv)
            if log
                println("          $k = ", current[k])
            end
        end
    end 
    
    if isnothing(constructor)
        return (; current...)
    else
        if usePath
            obj = Core.eval(modelModule, :($constructor(; path = $path, $current...))) 
        else
            obj = Core.eval(modelModule, :($constructor(; $current...)))
        end
        if log
            println("    +++ $path: typeof(obj) = ", typeof(obj), ", obj = ", obj, "\n\n")    
        end
        return obj        
    end
end




"""
    getIdParameter(par, id)
    
Search recursively in NamedTuple `par` for a NamedTuple that has 
`key = :_id, value = id` and return this NamedTuple or nothing,
if not present.
"""
function getIdParameter(par::NamedTuple, id::Int)
    if haskey(par, :_id) && par[:_id] == id
        return par
    else
        for (key,value) in zip(keys(par), par)
            if typeof(value) <: NamedTuple
                result = getIdParameter(value, id)
                if !isnothing(result)
                    return result
                end
            end
        end   
    end
    return nothing
end
