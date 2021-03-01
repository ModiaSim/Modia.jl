module TestTwoInertiasAndIdealGearWithUnits

using TinyModia
using DifferentialEquations
using ModiaPlot
using Unitful
using ModiaBase

TwoInertiasAndIdealGearWithUnits = Model(
    J1 = 0.0025u"kg*m^2",
    J2 = 170u"kg*m^2",
    r  = 105,
    tau_max = 1u"N*m",

    start = Map(phi2 = 0.5u"rad", w2 = 0.0u"rad/s", tau2 = 0u"N*m"), 
    
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

twoInertiasAndIdealGearWithUnits = @instantiateModel(TwoInertiasAndIdealGearWithUnits)


simulate!(twoInertiasAndIdealGearWithUnits, Tsit5(), stopTime = 4.0, 
          requiredFinalStates=[1.5628074713622309, -6.878080753044174e-5])
          
plot(twoInertiasAndIdealGearWithUnits, ["phi2", "w2", "der(w2)"])

end