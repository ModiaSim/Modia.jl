module TestFirstOrder

using TinyModia
using DifferentialEquations
using ModiaPlot

# using RuntimeGeneratedFunctions
# RuntimeGeneratedFunctions.init(@__MODULE__)

inputSignal(t) = sin(t)

FirstOrder = Model(
    T = 0.2,
    init = Map(x=0.3),
    equations = :[u = inputSignal(time/u"s"),
                  T * der(x) + x = u,
                  y = 2*x]
)

firstOrder = @instantiateModel(FirstOrder, logCode=true)

simulate!(firstOrder, Tsit5(), stopTime = 10, log=false, requiredFinalStates = [-0.3617373025974107])

plot(firstOrder, ["u", "x", "der(x)", "y"])

end