module TestFirstOrder

using TinyModia
using DifferentialEquations
using ModiaPlot

FirstOrder = Model(
    T = 0.2,
    u = input | Map(start=0),
    y = output,
    x = Var(init=0.3),
    equations = :[
        T * der(x) + x = u,
        y = 2*x]
)

firstOrder = @instantiateModel(FirstOrder, logModel=true, logDetails=true, log=true, logCode=false)

#simulate!(firstOrder, Tsit5(), stopTime = 10, log=false, requiredFinalStates = [-0.3617373025974107])

#plot(firstOrder, ["u", "x", "der(x)", "y"])

end