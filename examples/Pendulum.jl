module PendulumSimulation

using Modia, ModiaPlot, Measurements

Pendulum = Model(
   L = 0.8u"m",
   m = 1.0u"kg",
   d = 0.5u"N*m*s/rad",
   g = 9.81u"m/s^2",
   phi = Var(init = 1.57*u"rad"),
   w   = Var(init = 0u"rad/s"),
   equations = :[
          w = der(phi)
        0.0 = m*L^2*der(w) + d*w + m*g*L*sin(phi)
          r = [L*cos(phi), -L*sin(phi)]
   ]
)


pendulum1 = @instantiateModel(Pendulum)
simulate!(pendulum1, Tsit5(), stopTime = 10.0u"s", log=true, 
          requiredFinalStates = [-0.04808768814816799, -0.02484926198404018])
plot(pendulum1, [("phi", "w"); "r"], figure = 1)


PendulumWithUncertainties = Pendulum | Map(L = (0.8 ± 0.2)u"m",
                                           m = (1.0 ± 0.2)u"kg",
                                           d = (0.5 ± 0.2)u"N*m*s/rad")

pendulum2 =  @instantiateModel(PendulumWithUncertainties,
                               FloatType = Measurement{Float64})

simulate!(pendulum2, Tsit5(), stopTime = 10.0u"s", 
          requiredFinalStates = Measurements.Measurement{Float64}[-0.0480877836976501 ± 0.5039394725344267, -0.024849277709984924 ± 0.049301962043600225] )
plot(pendulum2, [("phi", "w"); "r"], figure = 2)


end