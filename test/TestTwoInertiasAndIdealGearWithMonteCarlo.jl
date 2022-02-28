module TestTwoInertiasAndIdealGearWithMonteCarlo

using Modia
@usingModiaPlot
using Modia.MonteCarloMeasurements
using Modia.MonteCarloMeasurements.Distributions

# The number of particles must be the same as for FloatType
const nparticles = 100
uniform(vmin,vmax) = StaticParticles(nparticles,Distributions.Uniform(vmin,vmax))

TwoInertiasAndIdealGearWithMonteCarlo = Model(
    J1 = 0.0025,
    J2 = uniform(50, 170),
    r  = 105,
    tau_max = 1,
    phi2 = Var(init = uniform(0.5,0.55)), 
    w2 = Var(init = 0.0),
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

twoInertiasAndIdealGearWithMonteCarlo = @instantiateModel(TwoInertiasAndIdealGearWithMonteCarlo, unitless=true,
                                                          FloatType = StaticParticles{Float64,nparticles})

simulate!(twoInertiasAndIdealGearWithMonteCarlo, Tsit5(), stopTime = 4.0,
          logParameters=true, logStates=true)
          
plot(twoInertiasAndIdealGearWithMonteCarlo, ["phi2", "w2"],
     MonteCarloAsArea=false)

end