#=
Propagates values in NamedTuples and instantiates dependent objects

* Developer: Hilding Elmqvist, Mogram AB (propagate)
*            Martin Otter, DLR (instantiateDependentObjects)
* First version: March 2021
* License: MIT (expat)

=#


using DataStructures: OrderedDict


subst(ex, environment) = ex
subst(ex::Symbol, environment) = if length(environment) == 0; ex elseif ex in keys(environment[end]); environment[end][ex] else subst(ex, environment[1:end-1]) end
subst(ex::Expr, environment) = Expr(ex.head, [subst(a, environment) for a in ex.args]...)

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



"""
    getIdParameter(par, id)
    
Search recursively in NamedTuple `par` for a NamedTuple that has 
`key = :_id, value = id` and return this NamedTuple or nothing,
if not present.
"""
function getIdParameter(par, id)
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


appendKey(path, key) = path == "" ? string(key) : path * "." * string(key)


"""
    instantiateDependentObjects(model)
    
Recursively traverse the NamedTuple `model` and instantiate all objects with
keys `_constructor = XXX` or `_constructor = (value = XXX, ...)` using
constructor `XXX`. Return the NamedTuple, where every NamedTuple with a `_constructor`
is replaced by the instantiated object (and `class, _constructor` removed in the
returned NamedTuple).
"""
function instantiateDependentObjects(modelModule, model, environment=[], path="")
    #println("\ninstantiateTypes of $path:", model)
    current = OrderedDict()
    constructor = nothing
    for (k,v) in zip(keys(model), model)
        #println("... for (k,v) in ...: $k = $v")
        if k == :_constructor
            if typeof(v) <: NamedTuple
                constructor = v[:value]
            else
                constructor = v
            end                
        elseif k == :class
            nothing;            
        elseif typeof(v) <: NamedTuple    
            current[k] = instantiateDependentObjects(modelModule, v, vcat(environment, [current]), appendKey(path, k))     
        elseif typeof(v) == Symbol
            current[k] = subst(v, vcat(environment, [current]))
        elseif Base.Meta.isexpr(v, :.) # v1.v2     
            v1 = v.args[1]
            v2 = v.args[2].value
            ev1 = subst(v1, vcat(environment, [current]))
            ev2 = subst(v2, [ev1])
            if typeof(ev2) <: NamedTuple
                current[k] = v
            else
                current[k] = ev2
            end
        else       
            current[k] = v
        end
    end 
    if isnothing(constructor)
        return (; current...)
    else
        #println("... constructor = $constructor, path = $path")
        return Core.eval(modelModule, :($constructor(; path = $path, $current...))) 
    end
end

