module TestFirstOrder

using Modia
@usingModiaPlot

inputSignal(t) = sin(t)

FirstOrder = Model(
    T = 0.2,
    x = Var(init=0.3),
    equations = :[u = inputSignal(time/u"s"),
                  T * der(x) + x = u,
                  y = 2*x]
)

firstOrder = @instantiateModel(FirstOrder, logCode=false)

simulate!(firstOrder, Tsit5(), stopTime = 10, log=true, requiredFinalStates = [-0.3617373025974107])

plot(firstOrder, ["u", "x", "der(x)", "y", "T"])

end