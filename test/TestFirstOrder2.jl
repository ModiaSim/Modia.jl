module TestFirstOrder2

using ModiaLang
using ModiaLang.DifferentialEquations
@usingModiaPlot
using Test

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

firstOrder = @instantiateModel(FirstOrder2, logCode=false)

simulate!(firstOrder, Tsit5(), stopTime = 10, merge = Map(T = 0.4, x = 0.9), 
          log=false, logParameters=true, logStates=true, 
          requiredFinalStates = [-0.17964872595554535])

# Test get_result(instantiatedModel)
println()
result1 = get_result(firstOrder)
@show(result1[1:10,:])
println()
@show(result1[1:10, ["time", "u", "y"]])

println()
result2 = get_result(firstOrder, onlyStates=true, extraNames=["y"])
@show(result2[1:10,:])

println()
result3 = get_result(firstOrder, extraNames=["y"])
@show(result3[1:10,:])

# Linearize
println("\n... Linearize at stopTime = 0 and 10:")
(A_0 , x_0)  = linearize!(firstOrder, analytic = true)
(A_10, x_10) = linearize!(firstOrder, stopTime=10, analytic = true) 
(A_10_numeric, x_10_numeric) = linearize!(firstOrder, stopTime=10, analytic=false) 
xNames = get_xNames(firstOrder)
@show xNames
@show A_0 , x_0
@show A_10, x_10
@show A_10_numeric, x_10_numeric
@test isapprox(A_0,[-1/0.4])
@test isapprox(A_0, A_10)


plot(result1, [("u", "x"), "der(x)", "y"])


FirstOrder3 = Model(
    T = 2u"hr",
    x = Var(init=1.0),
    equations = :[u = if after(1.5u"hr"); 1.0 else 0.0 end,
                  T * der(x) + x = u]
)
firstOrder3 = @instantiateModel(FirstOrder3, logCode=false)
simulate!(firstOrder3, Tsit5(), stopTime = 10u"hr")
plot(firstOrder3, [("u", "x"), "der(x)"], figure=2)

end