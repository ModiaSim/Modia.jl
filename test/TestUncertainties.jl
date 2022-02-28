module TestUncertainties

using Modia
@usingModiaPlot
using Modia.Measurements


TestModel = Model(
    T = 0.2 ± 0.02, # Cannot convert a particle distribution to a float if not all particles are the same.
    u0 = 10 ± 1, # Cannot convert a particle distribution to a float if not all particles are the same.
    k = 1 ± 0.1,
    x = Var(init = 5 ± 0.5),
    v = Var(init = 1 ± 0.1),

    m = 10 ± 1,
    F = 20 ± 2,
    equations = :[
    T*der(x) = u0
    y = k*x

    F = m * a
    der(v) = a # Cannot convert a particle distribution to a float if not all particles are the same.
    ]
)

model = @instantiateModel(TestModel, FloatType = Measurement{Float64})

simulate!(model, stopTime = 1.0, log=false,
          requiredFinalStates = Measurements.Measurement{Float64}[54.99999999999998 ± 7.08872343937891, 2.999999999999999 ± 0.3])

plot(model, ["T", "x", "der(x)", "y", "a", "der(v)"], figure=1)




include("../models/Electric.jl")

FilterCircuit = Model(
    R = Resistor        | Map(R = (100±10)u"Ω"),
    C = Capacitor       | Map(C = (0.01±0.001)u"F", v=Var(init=(0.0±1.0)u"V")),
    V = ConstantVoltage | Map(V = (10.0±1.0)u"V"),
    ground = Ground,
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n, ground.p)
    ]
)

filterCircuit = @instantiateModel(FilterCircuit, FloatType=Measurements.Measurement{Float64})
simulate!(filterCircuit, Tsit5(), stopTime = 3.0, logParameters = true, logStates = true) 
plot(filterCircuit, ("V.v", "C.v"), figure=2)          
   
        
Pendulum = Model(
   L = (0.8±0.1)u"m",
   m = (1.0±0.1)u"kg",
   d = (0.5±0.05)u"N*m*s/rad",
   g = 9.81u"m/s^2",
   phi = Var(init = (pi/2±0.1)*u"rad"), 
   w   = Var(init = 0u"rad/s"),
   equations = :[
          w = der(phi)
        0.0 = m*L^2*der(w) + d*w + m*g*L*sin(phi)
   ]
)

pendulum = @instantiateModel(Pendulum, FloatType=Measurements.Measurement{Float64})
simulate!(pendulum, Tsit5(), stopTime = 10.0)
plot(pendulum, "phi", figure=4)
        
          
end