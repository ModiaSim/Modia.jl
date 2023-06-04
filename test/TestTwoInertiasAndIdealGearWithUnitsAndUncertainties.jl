module TestTwoInertiasAndIdealGearWithUnitsAndUncertainties

using Modia
@usingModiaPlot
using Modia.Measurements


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

simulate!(twoInertiasAndIdealGear, Tsit5(), stopTime = 4.0, logParameters=true, logStates=true)

plot(twoInertiasAndIdealGear, ["phi2", "w2", "der(w2)"])

# Linearize
println("\n... Analytic linearization")
(A1, x1) = linearize!(twoInertiasAndIdealGear, stopTime=4, analytic = true)
xNames = get_xNames(twoInertiasAndIdealGear)
@show xNames
println(IOContext(stdout, :error_digits=>15), "A1 = ", A1, ", x1 = ", x1)

println("\n... Numeric linearization with Float64")
(A2, x2) = linearize!(twoInertiasAndIdealGear, stopTime=4, analytic=false)
println(IOContext(stdout, :error_digits=>15), "A2 = ", A2, ", x2 = ", x2)

#= DoubleFloats not defined
using Modia.Test
println("\n... Numeric linearization with Double64")
using DoubleFloats
twoInertiasAndIdealGear2 = InstantiatedModel{Measurement{Double64}}(twoInertiasAndIdealGear)
(A3, x3) = linearize!(twoInertiasAndIdealGear2, stopTime=3, analytic=false)
println(IOContext(stdout, :error_digits=>15), "A3 = ", A3, ", x3 = ", x3)
=#


end