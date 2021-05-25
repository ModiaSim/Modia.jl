module TestFirstOrder2

using TinyModia
using DifferentialEquations
using ModiaPlot

# using RuntimeGeneratedFunctions
# RuntimeGeneratedFunctions.init(@__MODULE__)

inputSignal(t) = sin(t)

FirstOrder1 = Model(
    T = 0.2,
    x = Var(init=0.3),
    equations = :[u = inputSignal(time/u"s"),
                  T * der(x) + x = u,
                  y = 2*x]
)

FirstOrder2 = FirstOrder1 | Map(T = 0.3, x = Var(init=0.6))

firstOrder = @instantiateModel(FirstOrder2, logCode=true)

simulate!(firstOrder, Tsit5(), stopTime = 10, merge = Map(T = 0.4, x = 0.9), 
          log=true, logParameters=true, logStates=true, 
          requiredFinalStates = [-0.17964872595554535])

# Test get_result(instantiatedModel)
println()
result1 = get_result(firstOrder)
@show(result1[1:10,:])

println()
result2 = get_result(firstOrder, onlyStates=true, extraNames=["y"])
@show(result2[1:10,:])

println()
result3 = get_result(firstOrder, extraNames=["y"])
@show(result3[1:10,:])

plot(result1, [("u", "x"), "der(x)", "y"])

end