module TestDeclarations

using Modia

struct Var
	var::NamedTuple
	Var(; kwargs...) = new((;kwargs...) )
end

struct Model
	model::NamedTuple
	Model(; kwargs...) = new((;kwargs...) )
end

struct Map
	map::NamedTuple
	Map(; kwargs...) = new((;kwargs...) )
end

Pin = Model( v = Var(potential=true, nominal=10), i = Var(;flow=true) )

@show Pin

Input(; kwargs...) = Var(;input=true, kwargs...)
Output(; kwargs...) = Var(;output=true, kwargs...)
Potential(; kwargs...) = Var(;potential=true, kwargs...)
Flow(; kwargs...) = Var(;flow=true, kwargs...)

using Base.Meta: isexpr
using Modia.OrderedCollections: OrderedDict
using ModiaBase.Symbolic: removeBlock


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

global logMerge = true 

function setLogMerge(val)
    global logMerge
    logMerge = val
end


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
    return (; gen...)
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

mergeModels(m1::Model, m2::Model, env=Symbol()) = recursiveMerge(m1.model, m2.model)
mergeModels(m1::Model, m2::Map, env=Symbol()) = recursiveMerge(m1.model, m2.map)

#Base.:|(m::Model, n::Model) =  recursiveMerge(m.model, n.model)
#Base.:|(m::Model, n::Map) = recursiveMerge(m.model, n.map)
#Base.:|(m, n::Map) = recursiveMerge(m, n.map)
Base.:|(m::Model, n::Model) =  (:MergeModel, m.model, n.model)
Base.:|(m::Model, n::Map) =  (:MergeMap, m.model, n.map)
Base.:|(m, n::Map) =  (:MergeMap, m, n.map)

Redeclare = ( _redeclare = true, )
Replace(i, e) = (:replace, i, e)
Final(x) = x

# ----------------------------------------------------------

@time Pin = Model( v = Potential(), i = Flow() )

OnePort = Model( p = Pin, n = Pin, equations = :[
        v = p.v - n.v
        0 = p.i + n.i
        i = p.i ] 
    )

@showModel(OnePort.model)

@time Resistor = OnePort | Model( R = Final(1.0u"Ω"), equations = Replace(1, :( R*i = v ) ))

@showModel(Resistor)

Capacitor = OnePort | Model( C = 1.0u"F", v = Map(init=0.0u"V"), equations = :[ C*der(v) = i ] ) 

Inductor = OnePort | Model( L = 1.0u"H", i = Map(init=0.0u"A"), equations = :[ L*der(i) = v ] )

ConstantVoltage = OnePort | Model( V = 1.0u"V", equations = :[ v = V ] )

Ground = Model( p = Pin, equations = :[ p.v = 0.0u"V" ] ) 

@time Filter = Model(
    R = Resistor | Map(R=0.5u"Ω"),
    C = Capacitor | Map(C=2.0u"F", v=Map(init=0.1u"V")),
    V = ConstantVoltage | Map(V=10.0u"V"),
    ground = Ground,
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n, ground.p)
    ]
)

@showModel(Filter.model)

end
