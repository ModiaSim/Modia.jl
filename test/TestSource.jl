module TestSource

using ModiaLang
using ModiaLang.DifferentialEquations
@usingModiaPlot



SineSource = Model(
    f = 2.0,
    equations = :[y = sin(2*3.14*f*time)]
)

sineSource = @instantiateModel(SineSource, unitless=true)

simulate!(sineSource, Tsit5(), stopTime = 1.0, log=false, logStates=false)

plot(sineSource, "y")

end