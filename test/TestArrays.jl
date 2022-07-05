module TestArrays

using Modia
using Modia.StaticArrays

@usingModiaPlot


LinearODE = Model(
    A = parameter | [-1.0  0.0;
                      0.0 -2.0],
    x = Var(init=[1.0, 2.0]),
    equations = :[der(x) = A*x]
)

linearODE = @instantiateModel(LinearODE, logCode=true)
simulate!(linearODE, stopTime = 2, log=false, logStates=false, requiredFinalStates=[0.13533533463680386, 0.036632273646086545])
plot(linearODE, ["x", "der(x)", "A"], figure=1, heading="LinearODE with length(x)=2")

simulate!(linearODE, stopTime = 2, log=false, logStates=false, requiredFinalStates=[0.14886882881582736, 0.038462894626776434, 0.00768439894426358], 
            merge = Map(A = [-1.0  0.0  0.0;
                              0.0 -2.0   0.0;
                              0.0  0.0 -3.0], x=[1.1, 2.1, 3.1]))
plot(linearODE, ["x", "der(x)"], figure=2, heading="LinearODE with length(x)=3")


LinearODE2 = Model(
    A = parameter | SMatrix{2,2}([-1.0  0.0;
                                   0.0 -2.0]),
    x = Var(init = SVector{2}(1.0, 2.0)),
    equations = :[der(x) = A*x]
)

linearODE2 = @instantiateModel(LinearODE2, logCode=true)
simulate!(linearODE2, stopTime = 2, log=false, logStates=false, requiredFinalStates=[0.13533533463680386, 0.036632273646086545])
@show typeof(linearODE2.evaluatedParameters[:A])
plot(linearODE2, ["x", "der(x)"], figure=3, heading="LinearODE2 with static arrays.")


end