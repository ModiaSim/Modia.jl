module TestUnitsAndUncertainites

using Modia
@usingModiaPlot
using Modia.Measurements

TestModel = Model(
    T = (0.2 ± 0.02)u"s",
    u0 = 8.5u"kg" ± 0.85u"kg",
    k = (1 ± 0.1)u"kg^-1",
    m = (10 ± 1)u"kg",
    F = (20 ± 2)u"N",

    x = Var(init = (5.5 ± 0.55)u"kg"),
    v = Var(init = (1 ± 0.1)u"m/s"), 
    r = Var(init = 1u"m/s"), 
    equations = :[
    T*der(x) + x = u0
    y = x
    z = x*k

    F = m * a
    der(v) = a
    der(r) = v ]
)

model = @instantiateModel(TestModel, FloatType = Measurement{Float64})

simulate!(model, Tsit5(), stopTime = 1.0)

plot(model, ["T", "x", "der(x)", "y", "a", "v", "r"])

end