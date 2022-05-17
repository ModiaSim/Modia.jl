module TestParameter

using Modia
@usingModiaPlot
using Modia.Test


include("$(Modia.path)/models/Blocks.jl")


inputSignal(t) = sin(t)

FirstOrder = Model(
    T = parameter | Var(value=0.2),
    x = Var(init=0.3),
    equations = :[u = inputSignal(time/u"s"),
                  T * der(x) + x = u,
                  y = 2*x]
)

firstOrder = @instantiateModel(FirstOrder)

simulate!(firstOrder, Tsit5(), stopTime = 10, log=false, logParameters=false, logEvaluatedParameters=false, requiredFinalStates = [-0.3617373025974107])

plot(firstOrder, ["u", "x", "der(x)", "y"], figure=1)



# Second order system:
#     w = 1.0
#     D = 0.1
#     k = 2.0
#     der(x1) = x2
#     der(x2) = -w^2*x1 - 2*D*w*x2 + w^2*u
#          y  = k*x1
SecondOrder = Model(
    w = 20.0,
    D =  0.1,
    k =  2.0,
    sys = StateSpace | Map(A = parameter | :([  0        1;
                                              -w^2  -2*D*w]),
                           B = parameter | :([0; w^2]),
                           C = parameter | :([k 0]),
                           D = parameter | :(zeros(1)),
                           x = Var(init = zeros(2))),
    equations = :[sys.u = 1.0]
)

@test_skip begin
secondOrder = @instantiateModel(SecondOrder, unitless=true)

simulate!(secondOrder, merge=Map(D=0.3), logParameters=true, logEvaluatedParameters=true,
          requiredFinalStates = [0.9974089572681231, 0.011808652820321723])
plot(secondOrder, [("sys.u", "sys.y"); "sys.x"], figure = 2)

simulate!(secondOrder, merge=Map(w=100.0, D=0.1, k=1.1), logParameters=true, logEvaluatedParameters=true,
          requiredFinalStates = [0.9999806306477742, -0.00391695258236727])
plot(secondOrder, [("sys.u", "sys.y"); "sys.x"], figure = 3)
end
end