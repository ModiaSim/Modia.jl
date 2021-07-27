module SimpleFilters

using Modia
@usingModiaPlot

SimpleModel = Model(
    T = 0.2,
    x = Var(init=0.5),
    equation = :[T * der(x) + x = 2]
)

# @showModel(SimpleModel)

model = @instantiateModel(SimpleModel)
simulate!(model, Tsit5(), stopTime = 5, requiredFinalStates = [1.9999996962023168])
plot(model, ["x"])

LowPassFilter = Model(
    T = 0.2,
    u = input,
    y = output | Var(:x),
    x = Var(init=0.0u"V"),
    equations = :[T * der(x) + x = u],
)

# @showModel(LowPassFilter)

HighPassFilter = LowPassFilter | Model(
        y = Var(:(-x + u)),
    )

# @showModel(HighPassFilter)

# setLogMerge(true)

LowAndHighPassFilter = LowPassFilter | Model(
    y = nothing,
    low = output | Var(:x),
    high = output | Var(:(-x + u)),
)

setLogMerge(false)

# @showModel(LowAndHighPassFilter)

TestLowAndHighPassFilter = LowAndHighPassFilter | Model(
    u = :(sin( (time+1u"s")*u"1/s/s" * time)*u"V"),
    x = Var(init=0.2u"V")
)

# @showModel(TestLowAndHighPassFilter)

model = @instantiateModel(TestLowAndHighPassFilter, log=false, logCode=false)
simulate!(model, Tsit5(), stopTime = 5, requiredFinalStates = [-0.22633046061014014])
plot(model, ["low", "high"])

TwoFilters = Model(
    high = HighPassFilter,
    low = LowPassFilter,
)

# @showModel(TwoFilters)

BandPassFilter = Model(
    u = input,
    y = output,
    high = HighPassFilter | Map(T=0.5, x=Var(init=0.1u"V")),
    low = LowPassFilter | Map(x=Var(init=0.2u"V")),
    equations = :[
        high.u = u,
        low.u = high.y,
        y = low.y]
)

# @showModel(BandPassFilter)

TestBandPassFilter = BandPassFilter | Model(
    u = :(sin( (0.05*time+1u"s")*u"1/s/s" * time)*u"V"),
)

# @showModel(TestBandPassFilter)

setLogMerge(false)

model = @instantiateModel(TestBandPassFilter, logDetails=false)
simulate!(model, Tsit5(), stopTime = 50, requiredFinalStates = [-0.25968254453053435, -0.6065869606784539])
plot(model, ["u", "y"])

end