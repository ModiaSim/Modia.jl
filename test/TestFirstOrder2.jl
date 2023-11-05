module TestFirstOrder2

using Modia
@usingModiaPlot
using Modia.Test

inputSignal(t) = sin(t)

FirstOrder1 = Model(
    T = 0.2u"s",
    x = Var(init=0.3),
    equations = :[u = inputSignal(time/u"s"),
                  T * der(x) + x = u,
                  y = 2*x]
)

FirstOrder2 = FirstOrder1 | Map(T = 0.3u"s", x = Var(init=0.6))

firstOrder = @instantiateModel(FirstOrder2, logCode=false)

simulate!(firstOrder, Tsit5(), stopTime = 10, merge = Map(T = 0.4u"s", x = 0.9), 
          log=false, logParameters=true, logStates=true, 
          requiredFinalStates = [-0.17964872595554535])

# Get result info
println("\n\n+++ Use SignalTables functions for post processing")

println("\n... Show overview of result")
showInfo(firstOrder)

println("\n... Get signal names, signals and signal info")
@show getSignalNames(firstOrder)
@show getStateNames(firstOrder)
sig_x = getSignal(firstOrder, "x")
@show getSignalInfo(firstOrder, "x")

println("\n... Get values")
@show getValues(firstOrder, "time")[1:5]
@show getValues(firstOrder, "y")[1:5]
@show getValue( firstOrder, "T")

@show getValuesWithUnit(firstOrder, "time")[1:5]
@show getValuesWithUnit(firstOrder, "y")[1:5]
@show getValueWithUnit( firstOrder, "T")

sig_der_x_flattened = getFlattenedSignal(firstOrder, "der(x)")
@show sig_der_x_flattened[:flattenedValues][1:5]
@show sig_der_x_flattened[:legend]

println("\n... Store result on file in JSON format")
writeSignalTable("TestFirstOrder2.json", firstOrder; indent=2, log=true)

println("\n... Store states on file in JSON format")
stateNames = getStateNames(firstOrder)
writeSignalTable("TestFirstOrder2_states.json", firstOrder; signalNames=stateNames, indent=2, log=true)


println("\n\n+++ Check deprecated functions\n")
result1 = get_result(firstOrder)
@show(result1[1:5,:])
println()
@show(result1[1:5, ["time", "u", "y"]])

println()
result2 = get_result(firstOrder, onlyStates=true, extraNames=["y"])
@show(result2[1:5,:])

println()
result3 = get_result(firstOrder, extraNames=["y"])
@show(result3[1:5,:])

@show signalNames(firstOrder)
@show timeSignalName(firstOrder)
@show hasOneTimeSignal(firstOrder)
printResultInfo(firstOrder)


# Linearize
println("\n\n+++ Linearize at stopTime = 0 and 10:")
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


plot(firstOrder, [("u", "x"), "der(x)", "y"])


FirstOrder3 = Model(
    T = 2u"hr",
    x = Var(init=1.0),
    equations = :[u = if after(1.5u"hr"); 1.0 else 0.0 end,
                  T * der(x) + x = u]
)
firstOrder3 = @instantiateModel(FirstOrder3, logCode=false)
#simulate!(firstOrder3, Tsit5(), stopTime = 10u"hr")
simulate!(firstOrder3, stopTime = 10u"hr")
plot(firstOrder3, [("u", "x"), "der(x)"], figure=2)

# Test all options
println("\n... Test all options of @instantiateModel(..)")
firstOrder3b = @instantiateModel(FirstOrder3, evaluateParameters=true, log=true, logModel=true, logDetails=true, logStateSelection=true,
                                logCode=true, logExecution=true, logCalculations=true, logTiming=true)


end