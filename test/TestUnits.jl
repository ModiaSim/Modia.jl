module TestUnits

using TinyModia
using DifferentialEquations
using ModiaPlot
using Test
using Unitful


UnitTest = Model(
    
    T = 0.2u"s",
    u0 = 8.5u"kg",
    k = 1u"kg^-1",
    m = 10u"kg",
    F = 20u"N",

    init = Map(x = 5500.0u"g",
        v = 1u"m/s"),
    equations = :[
        T*der(x) + x = u0
        y = x
        z = x*u"kg^-1"

        F = m * a
        der(v) = a ]
)

model = @instantiateModel(UnitTest)

simulate!(model, Tsit5(), stopTime = 1.0, requiredFinalStates = [5514.9625624219525, 2.9999999999999996])

plot(model, ["T", "x", "der(x)", "y", "a", "der(v)"])

end