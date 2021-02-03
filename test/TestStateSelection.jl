module TestStateSelection

using Unitful
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
    init = Map(phi=0.0u"rad"),
    equation = :[
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

drive = @instantiateModel(Drive, logDetails=true, logStateSelection=true)

end
