module TestInputOutput

using TinyModia
using DifferentialEquations
@usingModiaPlot

FirstOrder = Model(
    T = 0.2,
    u = input | Map(start=0),
    y = output,
    x = Var(init=0.3),
    equations = :[
        T * der(x) + x = u,
        y = 2*x]
)

firstOrder = @instantiateModel(FirstOrder, logModel=false, logDetails=false, log=true, logCode=true)

simulate!(firstOrder, Tsit5(), stopTime = 2, log=true, requiredFinalStates = [1.362130067500907e-5])

plot(firstOrder, ["u", "x", "der(x)", "y"])

end