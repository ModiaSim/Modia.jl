module TestTwoInertiasAndIdealGearWithUnitsAndMonteCarlo

using Modia
@usingModiaPlot
using Modia.MonteCarloMeasurements
using Modia.MonteCarloMeasurements.Distributions


# The number of particles must be the same as for FloatType
const nparticles = 100
uniform(vmin,vmax) = StaticParticles(nparticles,Distributions.Uniform(vmin,vmax))


TwoInertiasAndIdealGearWithUnitsAndMonteCarlo = Model(
    J1 = 0.0025u"kg*m^2",
    J2 = uniform(50.0, 170.0)u"kg*m^2",
    r  = 105.0,
    tau_max = 1.0u"N*m",
    
    phi2 = Var(start = 0.5u"rad"), 
    w2   = Var(start = 0.0u"rad/s"),
    tau2 = Var(start = 0u"N*m"),   

    equations = :[
        tau = if time < 1u"s"; tau_max elseif time < 2u"s"; 0.0u"N*m" elseif time < 3u"s"; -tau_max else 0.0u"N*m" end,

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

twoInertiasAndIdealGear = @instantiateModel(TwoInertiasAndIdealGearWithUnitsAndMonteCarlo,
                                            FloatType = StaticParticles{Float64,nparticles}, logCode=true)

simulate!(twoInertiasAndIdealGear, Tsit5(), stopTime = 4.0, logProgress=true)

plot(twoInertiasAndIdealGear, ["phi2", "w2"])

end