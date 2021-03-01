module TestTwoInertiasAndIdealGearWithMonteCarlo

using TinyModia
using DifferentialEquations
using ModiaPlot
using MonteCarloMeasurements
using Distributions
using ModiaBase
using Unitful

# The number of particles must be the same as for FloatType
const nparticles = 100
uniform(vmin,vmax) = StaticParticles(nparticles,Distributions.Uniform(vmin,vmax))

TwoInertiasAndIdealGearWithMonteCarlo = Model(
    J1 = 0.0025,
    J2 = uniform(50, 170),
    r  = 105,
    tau_max = 1,
    init = Map(phi2 = uniform(0.5,0.55), w2 = 0.0),
    equations = :[
        tau = if time < 1u"s"; tau_max elseif time < 2u"s"; 0 elseif time < 3u"s"; -tau_max else 0 end,
        
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

twoInertiasAndIdealGearWithMonteCarlo = @instantiateModel(TwoInertiasAndIdealGearWithMonteCarlo,
                                                          FloatType = StaticParticles{Float64,nparticles})

simulate!(twoInertiasAndIdealGearWithMonteCarlo, Tsit5(), stopTime = 4.0)
          
plot(twoInertiasAndIdealGearWithMonteCarlo, ["phi2", "w2"],
     MonteCarloAsArea=false)

end