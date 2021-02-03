module SimpleFilters

using TinyModia
using DifferentialEquations
using ModiaPlot
using Unitful

setLogMerge(false)

SimpleModel = Model(
    T = 0.2,
    equation = :(T * der(x) + x = 2),
)

# @showModel(SimpleModel)

LowPassFilter = Model(
    T = 0.2,
    inputs = :[u],
    outputs = :[y],
    init = Map(x=0),
    equations = :[T * der(x) + x = u],
    y = :x,
)

# @showModel(LowPassFilter)

HighPassFilter = merge(LowPassFilter, Model(
        y = :(-x + u),
    )
)

# @showModel(HighPassFilter)

LowAndHighPassFilter = LowPassFilter | Model(
    outputs = :[low, high],
    y = nothing,
    low = :x,
    high = :(-x + u),
)

# @showModel(LowAndHighPassFilter)

TestLowPassFilter = LowAndHighPassFilter | Model(
    w = :((time+1u"s")*u"1/s/s"),
    u = :(sin(w*time)*u"V"),
    init = Map(x=0.2u"V")
)

# @showModel(TestLowPassFilter)

TwoFilters = Model(
    high = HighPassFilter,
    low = LowPassFilter,
)

# @showModel(TwoFilters)

BandPassFilter = Model(
    inputs = :[u],
    outputs = :[y],
    high = HighPassFilter | Map(T=0.5, init=Map(x=0.1u"V")),
    low = LowPassFilter | Map(init=Map(x=0.2u"V")),
    equations = :[
        high.u = u,
        low.u = high.y,
        y = low.y]
)

# @showModel(BandPassFilter)

TestBandPassFilter = BandPassFilter | Model(
    w = :((0.05*time+1u"s")*u"1/s/s"),
    u = :(sin(w*time)*u"V"),
)

# @showModel(TestBandPassFilter)

setLogMerge(false)

model = @instantiateModel(TestBandPassFilter, logDetails=false)

simulate!(model, Tsit5(), stopTime = 50)
plot(model, ["u", "y"])

end