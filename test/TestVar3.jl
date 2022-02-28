module TestVar

using Modia
using Modia.Measurements
using Modia.StaticArrays

struct nnnVar
	var::NamedTuple
	Var(; kwargs...) = new((; kwargs...) )
end

struct Model
	model::NamedTuple
	Model(; kwargs...) = new((; kwargs...) )
end

struct Map
	map::NamedTuple
	Map(; kwargs...) = new((; kwargs...) )
end

Var(;kwargs...) = (;kwargs...)
Var(value; kwargs...) = (; value=value, kwargs...)

#Var(args...; kwargs...) = begin println("1"); (;args..., kwargs...) end
#Var(value, args...; kwargs...) = typeof(value) <: NamedTuple ? begin println("2", value, args, kwargs); (;value..., args..., kwargs...) end : begin println("3"); (;value = value, args..., kwargs...) end
xar(value, args...; kwargs...) = begin
	@show value args kwargs merge(value, args[1], args[2], (;kwargs...))
	end

var(value, args...; kwargs...) = typeof(value) <: NamedTuple ? 
	begin println("2", value, args, kwargs); value end : #merge(value, args, kwargs ) end : 
	begin println("3"); (;value = value, args..., kwargs...) end
#Var(value, args...; kwargs...) = (;value = value, args..., kwargs...)
#Var(value::Union{Number, String, Expr, AbstractArray}, args...; kwargs...) = (;value = value, args..., kwargs...)

#NotPair = Union{Number, String, Expr, AbstractArray} # ..., ....}
#Var(value::NotPair, args...; kwargs...) = (;value = value, args..., kwargs...)


Model(; kwargs...) = (; kwargs...)

constant = Var(constant = true)
parameter = Var(parameter = true)
input = Var(input = true)
output = Var(output = true)
potential = Var(potential = true)
flow = Var(flow = true)

range(min, max) = Var(min=min, max=max)

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
    return Var(; gen...)
end

Base.:|(m::Model, n::Model) =  recursiveMerge(m.model, n.model)
Base.:|(m::Model, n::Map) = recursiveMerge(m.model, n.map)
Base.:|(m, n::Map) = recursiveMerge(m, n.map)
Base.:|(m, n) = if !(typeof(n) <: NamedTuple); recursiveMerge(m, (; value=n)) else recursiveMerge(n, (value=m,)) end

#Base.:|(m::Model, n::Model) =  (:MergeModel, m.model, n.model)
#Base.:|(m::Model, n::Map) =  (:MergeMap, m.model, n.map)
#Base.:|(m, n::Map) =  (:MergeMap, m, n.map)
#Base.:|(m::Var, n::Var) =  (:MergeModel, m.var, n.var)

# ---------------------------------------------------------------------

#v = Var(potential, nominal=10)
#@show v

v1 = input | output | flow | Var(min=0, nominal=10)
@show v1

v2 = parameter | 5
@show v2

v3 = 10u"m/s" | parameter
@show v3

v4 = parameter | 5.5 | Var(min=0)
@show v4

v5 = parameter | Var(5, min=0)
@show v5

Pin = Model( v = potential | Var(nominal=10), i = flow )
@show Pin

v6 = Var(2.0)
@show v6

v7 = Var(2.0 ± 0.1)
@show v7

v8 = Var(SVector(1.0,2.0,3.0))
@show v8

v9 = v4 | Var(parameter=false)
@show v9

v10 = v4 | Var(parameter=nothing)
@show v10

v11 = v4 | range(0,10)
@show v11




end