module TestSource

using TinyModia
using DifferentialEquations
using ModiaPlot



SineSource = Model(
    f = 2.0,
    equations = :[y = sin(2*3.14*f*time)]
)

sineSource = @instantiateModel(SineSource, unitless=true)

simulate!(sineSource, Tsit5(), stopTime = 1.0, log=true, logStates=true)

plot(sineSource, "y")

end