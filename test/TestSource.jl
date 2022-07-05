module TestSource

using Modia
@usingModiaPlot



SineSource = Model(
    f = 2.0,
    equations = :[y = sin(2*3.14*f*time)]
)

sineSource = @instantiateModel(SineSource, unitless=true)

simulate!(sineSource, Tsit5(), stopTime = 1.0, log=false, logStates=true)

plot(sineSource, "y")

end