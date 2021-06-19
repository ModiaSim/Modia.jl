module Tutorial

using Modia
@usingModiaPlot

# 1 Getting started

using Modia
@usingModiaPlot

# Define model
SimpleModel = Model(
    T = 0.4,
    x = Var(init=0.2),
    equation = :[T * der(x) + x = 1],
)

# Transform to ODE form
simpleModel = @instantiateModel(SimpleModel)

# Simulate with a default integrator of DifferentialEquations
simulate!(simpleModel, stopTime = 1.2)

# Simulate with a specific integrator (Tsit5) and use a unit for stopTime
simulate!(simpleModel, Tsit5(), stopTime = 1.2u"s")

# Produce a line plot with GLMakie
plot(simpleModel, ("x", "der(x)"), figure=1)



# 2.1 - Equation oriented models

using Modia

LowPassFilter = Model(
    T = 0.2,
    u = input,
    y = output | Var(:x),
    x = Var(init=0),
    equation = :[T * der(x) + x = u],
)
    

   
# 2.2 Merging models  

HighPassFilter = LowPassFilter | Model( y = :(-x + u) )

setLogMerge(true)

LowAndHighPassFilter = LowPassFilter | Model(
        y = nothing,
        low = output | Var(:x),
        high = output | Var(:(-x + u)),
    )
    
setLogMerge(false)

@showModel LowAndHighPassFilter



# 2.3 Functions and tables

TestLowAndHighPassFilter = LowAndHighPassFilter | Model(
        u = :(sin( (time+1u"s")*u"1/s/s" * time)*u"V"),
        x = Var(init=0.2u"V")
    )
    
    
using Interpolations

table = CubicSplineInterpolation(0:0.5:2.0, [0.0, 0.7, 2.0, 1.8, 1.2])
TestLowAndHighPassFilter2 = TestLowAndHighPassFilter | Map(u = :(table(time*u"1/s")*u"V")) 
    
    
    
# 2.4 Hierarchical modeling

TwoFilters = (
    high = HighPassFilter,
    low = LowPassFilter,
)

BandPassFilter = (
    u = input,
    y = output,
    high = HighPassFilter | Map(T=0.5, x=Var(init=0.1u"V")),
    low = LowPassFilter | Map(x=Var(init=0.2u"V")),
    equations = :[
        high.u = u,
        low.u = high.y,
        y = low.y]
)

TestBandPassFilter = BandPassFilter | Map(
        u = :(sin( (time+1u"s")*u"1/s/s" * time)*u"V")
    )

#bandPassFilter = @instantiateModel(TestBandPassFilter, logStateSelection=true)
#simulate!(bandPassFilter, Tsit5(), stopTime = 50u"s")
#plot(bandPassFilter, ["u", "y"], figure=2)



end