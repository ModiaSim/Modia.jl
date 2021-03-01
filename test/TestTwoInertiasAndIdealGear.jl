module TestTwoInertiasAndIdealGear

using TinyModia
using DifferentialEquations
using ModiaPlot
using Unitful


TwoInertiasAndIdealGearTooManyInits = Model(
    J1 = 0.0025,
    J2 = 170,
    r  = 105,
    tau_max = 1,
    init = Map(phi1 = 0.0, w1 = 1.0, phi2 = 0.5, w2 = 0.0),
    equations = :[
        tau = if time < 1u"s"; tau_max elseif time < 2u"s"; 0 elseif time < 3u"s"; -tau_max else 0 end,

        # inertia1
        w1 = der(phi1),
        J1*der(w1) = tau - tau1,

        # gear
        phi1   = r*phi2,
        r*tau1 = tau2,

        # inertia2]
        w2 = der(phi2),
        J2*der(w2) = tau2
    ]
)

TwoInertiasAndIdealGear = TwoInertiasAndIdealGearTooManyInits | Map(init = Map(phi1=nothing, w1=nothing))

twoInertiasAndIdealGearTooManyInits = @instantiateModel(TwoInertiasAndIdealGearTooManyInits)
twoInertiasAndIdealGear             = @instantiateModel(TwoInertiasAndIdealGear)

println("Next simulate! should result in an error:\n")
simulate!(twoInertiasAndIdealGearTooManyInits, Tsit5(), stopTime = 4.0, log=true)

simulate!(twoInertiasAndIdealGear, Tsit5(), stopTime = 4.0, log=true,
          requiredFinalStates=[1.5628074713622309, -6.878080753044174e-5])
          
plot(twoInertiasAndIdealGear, ["phi2", "w2"])


end