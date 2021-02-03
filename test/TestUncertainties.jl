module TestUncertainties

using TinyModia
using DifferentialEquations
using ModiaPlot

using Measurements

TestModel = Model(
    T = 0.2 ± 0.02, # Cannot convert a particle distribution to a float if not all particles are the same.
    u0 = 10 ± 1, # Cannot convert a particle distribution to a float if not all particles are the same.
    k = 1 ± 0.1,
    init = Map(x = 5 ± 0.5,
        v = 1 ± 0.1),

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

plot(model, ["T", "x", "der(x)", "y", "a", "der(v)"])

end