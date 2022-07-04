module TestFirstOrder2

using Modia
@usingModiaPlot
using Modia.Test

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
showInfo(firstOrder)


# Get result info
println("\n... Get values")
t = getValues(firstOrder, "time")
y = getValues(firstOrder, "y")
T = getValue( firstOrder, "T")
tWithUnit = getValuesWithUnit(firstOrder, "time")
@show t[1:10]
@show y[1:10]
@show T
@show tWithUnit[1:10]

# Store result info on file
println("\n... Store result on file in JSON format")
writeSignalTable("TestFirstOrder2.json", firstOrder; indent=2, log=true)

println("\n... Store states on file in JSON format")
stateNames = getStateNames(firstOrder)
writeSignalTable("TestFirstOrder2_states.json", firstOrder; signalNames=stateNames, indent=2, log=true)

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


plot(firstOrder, [("u", "x"), "der(x)", "y"])


FirstOrder3 = Model(
    T = 2u"hr",
    x = Var(init=1.0),
    equations = :[u = if after(1.5u"hr"); 1.0 else 0.0 end,
                  T * der(x) + x = u]
)
firstOrder3 = @instantiateModel(FirstOrder3, logCode=false)
simulate!(firstOrder3, Tsit5(), stopTime = 10u"hr")
plot(firstOrder3, [("u", "x"), "der(x)"], figure=2)

# Test all options
println("\n... Test all options of @instantiateModel(..)")
firstOrder3b = @instantiateModel(FirstOrder3, evaluateParameters=true, log=true, logModel=true, logDetails=true, logStateSelection=true,
                                logCode=true, logExecution=true, logCalculations=true, logTiming=true)


end