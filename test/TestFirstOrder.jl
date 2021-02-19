module TestFirstOrder

using TinyModia
using DifferentialEquations
using ModiaPlot
using Unitful
using ModiaBase
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

FirstOrder = Model(
    T = 0.2,
    init = Map(x=0.3),
    equations = :[u = 1.0,
                  T * der(x) + x = u,
                  y = 2*x]
)

firstOrder = @instantiateModel(FirstOrder, logCode=true)

simulate!(firstOrder, Tsit5(), stopTime = 1.0, log=false,
          requiredFinalStates = [0.9952834203188597])

plot(firstOrder, ["T", "x", "der(x)", "y"])

end