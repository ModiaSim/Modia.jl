module TestUnits

using Modia
@usingModiaPlot


UnitTest = Model(
    
    T = 0.2u"s",
    u0 = 8.5u"kg",
    k = 1u"kg^-1",
    m = 10u"kg",
    F = 20u"N",

    x = Var(init = 5.5u"kg"),    # 5500.0u"g"),
    v = Var(init = 1u"m/s"),
    equations = :[
        T*der(x) + x = u0
        y = x
        z = x*u"kg^-1"

        F = m * a
        der(v) = a ]
)

model = @instantiateModel(UnitTest)

simulate!(model, Tsit5(), stopTime = 1.0, requiredFinalStates = [8.479786015016273, 2.999999999999999]) 

#plot(model, ["T", "x", "der(x)", "y", "a", "der(v)"])
plot(model, ["x", "der(x)", "y", "a", "der(v)"])

end