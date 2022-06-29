module TestStateSelection

using Modia
using Modia.Test

function checkStateSelection(model, x_names, linearEquations=[])
    # Check names of the states
    @show model.equationInfo
    
    for (i, state_i) in enumerate(model.equationInfo.x_info)
        @test state_i.x_name == x_names[i]
    end

    # Check linear equation systems
    if length(linearEquations) > 0
        @test model.equationInfo.linearEquations == linearEquations
    end
end


@testset "\nTest Modia/src/StateSelection.jl" begin

    @testset "... Test FirstOrder" begin
        FirstOrder = Model(
            T = 0.2,
            x = Var(init=0.3),
            equations = :[
                    u = 1.0
                der(x) = (u - x)/T
                    y = 2*x]
        )
        firstOrder = @instantiateModel(FirstOrder)
        checkStateSelection(firstOrder, ["x"])
    end
    

    @testset "... Test TwoCoupledInertias" begin
        TwoCoupledInertias = Model(
            J1_J = 1,
            J2_J = 1,
            J1_phi = Var(start=1.0), 
            J1_w   = Var(start=0.0), 
            J2_phi = Var(start=0.0), 
            J2_w   = Var(start=0.0),
            equations = :[
                J1_w = der(J1_phi)
                J1_J * der(J1_w) = J1_tau

                J2_w = der(J2_phi)
                J2_J * der(J2_w) = J2_tau

                J2_phi = J1_phi
                J2_tau + J1_tau = 0]
        )
        twoCoupledInertias = @instantiateModel(TwoCoupledInertias)
        checkStateSelection(twoCoupledInertias, ["J1_w", "J1_phi"], [(["der(J1_w)"], [], [1], 1, 1)])
    end


    @testset "... Test ODE with linear equations 1" begin
        ODEwithLinearEquations1 = Model(
            p1 = 4.0,
            p2 = 2,
            x6 = Var(start=1.0),
            equations = :[
                         x1 = sin(time)
                      p2*x2 = p1*x1
                2*x3 + 4*x4 = x1
                    x3 + x5 = 2*x2
                4*x4 + 5*x5 = x1 + x2
                    der(x6) = -x5*x6]
        )
        oDEwithLinearEquations1 = @instantiateModel(ODEwithLinearEquations1, unitless=true)
        checkStateSelection(oDEwithLinearEquations1, ["x6"], [(["x5"], [], [1], 1, 1)])
    end


    @testset "... Test ODE with linear equations 2" begin
        ODEwithLinearEquations2 = Model(
            p1=4,
            p2=2,
            x6 = Var(start=1.0),
            equations = :[
                        x1 = sin(time)
                     p2*x2 = p1*x1
                p2*x3 + p1*x4 + 5*der(x6) = x1
                    x3 + x5 = 2*x2
                4*x4 + 5*x5 = x1 + x2
                  3*der(x6) = -5*x3 - 2*x4]
        )
        oDEwithLinearEquations2 = @instantiateModel(ODEwithLinearEquations2, unitless=true)
        checkStateSelection(oDEwithLinearEquations2, ["x6"], [(["x5"], [], [1], 1, 1)])
    end


    @testset "... Test Multi-index DAE" begin
        MultiIndexDAE = Model(
            x2 = Var(init = 0.0),
            x2d = Var(init = 0.0),
            equations = :[
                u1 = sin(1*time)
                u2 = sin(2*time)
                u3 = sin(3*time)
                u4 = sin(4*time)
                u5 = sin(5*time)
                u6 = sin(6*time)
                u7 = sin(7*time)
                0 = u1 + x1 - x2
                0 = u2 + x1 + x2 - x3 + der(x6)
                0 = u3 + x1 + der(x3) - x4
                x1d  = der(x1)
                x2d  = der(x2)
                x3d  = der(x3)
                x6d  = der(x6)
                x6dd = der(x6d)
                0 = u4 + 2*der(x1d) + der(x2d) + der(x3d) + der(x4) + der(x6dd)
                0 = u5 + 3*der(x1d) + 2*der(x2d) + x5
                0 = u6 + 2*x6 + x7
                0 = u7 + 3*x6 + 4*x7]
        )
        multiIndexDAE = @instantiateModel(MultiIndexDAE, unitless=true)
        #@show multiIndexDAE.equationInfo.linearEquations
        checkStateSelection(multiIndexDAE, ["x2", "x2d"],
                            [(["x7"], [], [1], 1, 1),
                             (["der(x7)"], [], [1], 1, 1),
                             (["der(der(x7))"], [], [1], 1, 1),
                             (["der(der(der(x7)))"], [], [1], 1, 1),
                             (["der(x2d)"], [], [1], 1, 1)])
    end


    @testset "... Test free flying mass" begin
        FreeFlyingMass = Model(
            m =1.0,
            f=[1.0,2.0,3.0],
            r = Var(init=[0.1, 0.2, 0.3]), 
            v = Var(init=[-0.1, -0.2, -0.3]),
            equations = :[
                v = der(r)
                m*der(v) = f]
        )
        freeFlyingMass = @instantiateModel(FreeFlyingMass)
        checkStateSelection(freeFlyingMass, ["v", "r"], [])
    end

#= Does not yet work, because Pantelides does not yet support arrays 
    @testset "... Test sliding mass" begin
        SlidingMass = Model(
            m = 1.0,
            c = 1e4,
            d = 10.0,
            g = 9.81,
            n = [1.0, 1.0, 1.0],
            s = Var(init=1.0),
            equations = :[
                r = n*s
                v = der(r)
                m*der(v) = f + m*g*[0,-1,0] + u
                0 = n*f
                u = -(c*s + d*der(s))*n]
        )
        slidingMass = @instantiateModel(SlidingMass, logStateSelection=true)
        @show slidingMass.equationInfo.linearEquations
        #checkStateSelection(slidingMass, ["x6"], [(["x5"], [1], 1, 1)])
    end
=#
end


include("../models/AllModels.jl")

#=
Gear = Model(
    flange_a = Flange,
    flange_b = Flange,
    gear     = IdealGear_withSupport | Map(ratio = 105.0),
    fixed    = Fixed,
    spring   = Spring | Map(c=5.84e5u"N*m/rad"),
    damper1  = Damper | Map(d=500.0u"N*m*s/rad"),
    damper2  = Damper | Map(d=100.0u"N*m*s/rad"),

    connect = :[
        (flange_a     , gear.flange_a)
        (fixed.flange , gear.support, damper2.flange_b)
        (gear.flange_b, damper2.flange_a, spring.flange_a, damper1.flange_a)
        (flange_b     , spring.flange_b, damper1.flange_b)]
)

Drive = Model(
    inertia1 = Inertia,
    inertia2 = Inertia,
    gear     = Gear,
    equations = :[
        inertia1.flange_a.tau = 0u"N*m"
        inertia2.flange_b.tau = 0u"N*m"],
    connect = :[
        (inertia1.flange_b, gear.flange_a),
        (gear.flange_b    , inertia2.flange_a)
    ]
)
=#

# Linear 1D rotational damper
AbsoluteDamper = Model(
    flange = Flange,
    d = 1.0u"N*m*s/rad", # (info = "Damping constant"),
    phi = Var(init=0.0u"rad"),
    equations = :[
        phi = flange.phi
        w   = der(phi)
        flange.tau = d * w ]
)


Drive = Model(
    inertia = Inertia,
    damper  = AbsoluteDamper,
    equations = :[
        inertia.flange_a.tau = 0u"N*m"],
    connect = :[
        (inertia.flange_b, damper.flange)]
)
Drive2 = Drive | Map(damper = Map(phi=Var(init=1.0u"rad")))

drive1 = @instantiateModel(Drive)
drive2 = @instantiateModel(Drive2)

simulate!(drive1, Tsit5(), stopTime = 4.0)

println("Next simulate! should result in an error:\n")
simulate!(drive2, Tsit5(), stopTime = 4.0)


end
