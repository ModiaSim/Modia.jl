module TestTwoInertiasAndIdealGearWithUnits

using Modia
@usingModiaPlot
using Modia.Test

TwoInertiasAndIdealGearWithUnits = Model(
    J1 = 0.0025u"kg*m^2",
    J2 = 170u"kg*m^2",
    r  = 105.0,
    tau_max = 1u"N*m",

    phi2 = Var(start = 0.5u"rad"), 
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

twoInertiasAndIdealGearWithUnits = @instantiateModel(TwoInertiasAndIdealGearWithUnits, logCode=true)


simulate!(twoInertiasAndIdealGearWithUnits, Tsit5(), stopTime = 4.0, 
          logParameters=true, logStates=true,
          requiredFinalStates=[1.5628074713622309, -6.878080753044174e-5])
          
plot(twoInertiasAndIdealGearWithUnits, ["phi2", "w2", "der(w2)"])

# Linearize
println("\n... Linearize at stopTime = 0 and 4")
(A0, x0) = linearize!(twoInertiasAndIdealGearWithUnits, stopTime=0, analytic = true)
(A1, x1) = linearize!(twoInertiasAndIdealGearWithUnits, stopTime=4, analytic = true) 
xNames = get_xNames(twoInertiasAndIdealGearWithUnits)
@show xNames
@show A0, x0
@show A1, x1
@test isapprox(A0,[0.0 1.0; 0.0 0.0])
@test isapprox(A0, A1)

end