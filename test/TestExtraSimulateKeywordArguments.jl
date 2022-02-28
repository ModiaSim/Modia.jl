module TestExtraSimulateKeywordArguments

using Modia
@usingModiaPlot

# Register extra simulate! keyword arguments
registerExtraSimulateKeywordArguments([:dummy1, :dummy2, :dummy3])

function inputSignal(instantiatedModel, t)
    if isInitial(instantiatedModel)
        dict = get_extraSimulateKeywordArgumentsDict(instantiatedModel)
        dummy1 = get(dict, :dummy1, 1.0)
        dummy2 = get(dict, :dummy2, true)
        dummy3 = get(dict, :dummy3, "this is dummy3")
        println("\nExtra simulate! keyword arguments: dummy1 = $dummy1, dummy2 = $dummy2, dummy3 = \"$dummy3\"\n")
    end
    
    sin(t)
end

FirstOrder = Model(
    T = 0.2,
    x = Var(init=0.3),
    equations = :[u = inputSignal(instantiatedModel, time/u"s"),
                  T * der(x) + x = u,
                  y = 2*x]
)

firstOrder = @instantiateModel(FirstOrder, logCode=false)

simulate!(firstOrder, Tsit5(), stopTime = 10, log=false, dummy1 = 2.0, dummy3 = "value changed",
          requiredFinalStates = [-0.3617373025974107])

plot(firstOrder, ["u", "y"])

end