module TestMechanics

println("TestMechanics: Tests how 3D mechanics could be combined with Modia.")

using Modia
@usingModiaPlot

# Modia models
include("$(Modia.path)/models/Blocks.jl")
include("$(Modia.path)/models/Electric.jl")
include("$(Modia.path)/models/Rotational.jl")


# Generic Modia3D model stubs (always available)
RevoluteStub = Model(
    flange = Flange,

    equations = :[phi  = flange.phi 
                  w    = der(phi)
                 der_w = der(w)
                 tau   = flange.tau]
)


# Generic Modia3D definitions (always available)
mutable struct InstantiatedModel3D
    # Internal memory of 3D model
    initial::Bool

    InstantiatedModel3D() = new(true)
end


# Function ..._f1 for 3D mechanics
# function Pendulum_f1(_m::InstantiatedModel3D, phi::Float64, w::Float64, tau::Float64, 
#                       m::Float64, L::Float64, g::Float64)
function Pendulum_f1(_m::TestMechanics.InstantiatedModel3D, phi, w, tau, m, L, g)
#function Pendulum_f1(phi, w, tau, m, L, g)

#=
    if _m.initial
        println("... SimplePendulum_f1 called during initialization")
        _m.initial = false
    end
=#
    der_w = (tau - m*g*L*cos(phi))/m*L*L
    return der_w
end


# Modia model for the system
Pendulum = Model(
    model3D = InstantiatedModel3D(),
    
    # Parameters
    m = 1.0u"kg",
	L = 0.5u"m",
	g = 9.81u"m/s^2",

    # 3D model stubs
    rev = RevoluteStub | Map(init = Map(phi=1.5u"rad", w=2u"rad/s")),

    # Standard Modia models
    damper  = Damper | Map(d=0.4u"N*m*s/rad"),
	support = Fixed,

    connect = :[
        (rev.flange     , damper.flange_b),
        (damper.flange_a, support.flange)],
        
    equations = :(rev.der_w = Pendulum_f1(model3D, rev.phi, rev.w, rev.tau, m, L, g))
#    equations = :(rev.der_w = Pendulum_f1(rev.phi, rev.w, rev.tau, m, L, g))
)

pendulum = @instantiateModel(Pendulum, logCode=true, logExecution=false, unitless=true)

simulate!(pendulum, Tsit5(), stopTime=20.0, log=true)
plot(pendulum, ["rev.phi", "rev.w"])

end

