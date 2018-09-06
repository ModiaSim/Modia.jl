

using .Instantiation
using .ModelElaboration

using Unitful
using .StructuralTransform
using .Synchronous: sample, Clock, previous, hold, positive, positiveChange, positiveEdge

"""
Shortcut for `Variable`
"""
Var(; args...) = Variable(; args...)

"""
Create a floating-point `Variable`
"""
Float(value=nothing; info="", size=nothing, unit=NoUnits, displayUnit=NoUnits, 
    min=nothing, max=nothing, start=nothing, fixed::Bool=false, nominal=nothing,
    variability=continuous, 
    flow::Bool=false, state::Bool=true) = Variable(variability, Float64, size, value, 
    unit, displayUnit, min, max, start, fixed, nominal, info, flow, state, general)

"""
Create a boolean `Variable`
"""
Boolean(value=nothing; info="", size=nothing, unit=NoUnits, displayUnit=NoUnits, 
    min=nothing, max=nothing, start=nothing, fixed::Bool=false, nominal=nothing,
    variability=continuous, 
    flow::Bool=false, state::Bool=true) = Variable(variability, Bool, size, value, 
    unit, displayUnit, min, max, start, fixed, nominal, info, flow, state, general)

"""
Create an integer `Variable`
"""
Integ(value=nothing; info="", size=nothing, unit=NoUnits, displayUnit=NoUnits, 
    min=nothing, max=nothing, start=nothing, fixed::Bool=false, nominal=nothing,
    variability=continuous, 
    flow::Bool=false, state::Bool=true) = Variable(variability, Int, size, value, 
    unit, displayUnit, min, max, start, fixed, nominal, info, flow, state, general)

"""
Create a string `Variable`
"""
Str(value=nothing; info="", size=nothing, unit=NoUnits, displayUnit=NoUnits, 
    min=nothing, max=nothing, start=nothing, fixed::Bool=false, nominal=nothing,
    variability=continuous, 
    flow::Bool=false, state::Bool=true) = Variable(variability, String, size, value, 
    unit, displayUnit, min, max, start, fixed, nominal, info, flow, state, general)

#=
Float(; args...) = Var(T=Float64; args...)
Float0(; args...) = Float(size=(); args...)

Boolean(; args...) = Var(T=Bool; args...)
Boolean0(; args...) = Boolean(size=(); args...)

Integ(; args...) = Var(T=Int; args...)
Integ0(; args...) = Integ(size=(); args...)

Str(; args...) = Var(T=String; args...)
Str0(; args...) = Str(size=(); args...)
=#

"""
Create a `Variable` with `parameter` variability, meaning it 
is an input variable that is constant with time
"""
Parameter(; args...) = Variable(variability=parameter; args...)
Parameter(value; args...) = Variable(variability=parameter, value=value; args...)

"""
Shortcut for `Parameter`
"""
Par(; args...) = Variable(variability=parameter; args...)
Par(value; args...) = Variable(variability=parameter, value=value; args...)

"""
A value meant to be filled in later
"""
const undefined = nothing

