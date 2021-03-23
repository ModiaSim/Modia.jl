module TestTwoInertiasAndIdealGearWithUnitsAndUncertainties


using TinyModia
using DifferentialEquations
using ModiaPlot
using Measurements
using Unitful



TwoInertiasAndIdealGear = Model(
    J1 = 0.0025u"kg*m^2",
    J2 = (150.0 ± 20.0)u"kg*m^2",
    r  = 105.0,
    tau_max = 1u"N*m",

    phi2 = Var(start = (0.5 ± 0.05)u"rad"), 
    w2   = Var(start = 0.0u"rad/s"), 
    tau2 = Var(start = 0u"N*m"), 
    
    equations = :[
        tau = if time < 1u"s"; tau_max elseif time < 2u"s"; 0u"N*m" elseif time < 3u"s"; -tau_max else 0u"N*m" end,
    
        # inertia1
        w1 = der(phi1),
        J1*der(w1) = tau - tau1,
 
        # gear
        phi1   = r*phi2,
        r*tau1 = tau2,
    
        # inertia2
        w2 = der(phi2),
        J2*der(w2) = tau2
    ]
)

twoInertiasAndIdealGear = @instantiateModel(TwoInertiasAndIdealGear, FloatType = Measurement{Float64})

simulate!(twoInertiasAndIdealGear, Tsit5(), stopTime = 4.0)

plot(twoInertiasAndIdealGear, ["phi2", "w2", "der(w2)"])

end