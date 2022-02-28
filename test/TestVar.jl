module TestVar

using Modia.Measurements
using Modia.StaticArrays

Var(;kwargs...) = (;kwargs...)
Var(value; kwargs...) = (;value=value, kwargs...)

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

# ---------------------------------------------------------------------

#v = Var(potential, nominal=10)
#@show v

v = input | output | flow | Var(min=0, nominal=10)
@show v

v = parameter | Var(5, min=0)
@show v

v = parameter | Var(5, min=0)
@show v

Pin = Model( v = potential | Var(nominal=10), i = flow )
@show Pin

v2 = Var(2.0)
@show v2

v3 = Var(2.0 ± 0.1)
@show v3

v4 = Var(SVector(1.0,2.0,3.0))
@show v4

end