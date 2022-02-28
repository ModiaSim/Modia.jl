module TestStateSpace

using Modia
@usingModiaPlot

include("$(Modia.path)/models/Blocks.jl")


# Second order system:
#     w = 1.0
#     D = 0.1
#     k = 2.0
#     der(x1) = x2
#     der(x2) = -w^2*x1 - 2*D*w*x2 + w^2*u
#          y  = k*x1
SecondOrder1 = Model(
    w = 20.0,
    D =  0.1,
    k =  2.0,
    sys = StateSpace | Map(A = :([  0        1;
                                 -w^2  -2*D*w]),
                           B = :([0; w^2]),
                           C = :([k 0]),
                           D = :(zeros(1)),
                        x = Var(init = zeros(2)) ),
    equations = :[sys.u = 1.0]
)


secondOrder1 = @instantiateModel(SecondOrder1, unitless=true)
simulate!(secondOrder1, stopTime=1.0, merge=Map(D=0.3),
          requiredFinalStates = [0.9974089572681231, 0.011808652820321723])
plot(secondOrder1, [("sys.u", "sys.y"); "sys.x"], figure = 1)



SecondOrder2 = Model(
    sys = StateSpace | Map(   A = [   0.0    1.0;
                                   -100.0  -12.0],
                              B = [0; 400.0],
                              C = [2.0 0],
                              D = zeros(1,1),
                           x = Var(init = zeros(2))),
    equations = :[sys.u = 1.0]
)

secondOrder2 = @instantiateModel(SecondOrder2, unitless=true)
simulate!(secondOrder2, Tsit5(), stopTime=1.0, merge=Map(sys=Map(A = [0.0 1.0; -400.0  -12.0])),
          requiredFinalStates = [0.9974089572681231, 0.011808652820321723])
plot(secondOrder2, [("sys.u", "sys.y"); "sys.x"], figure = 2)

simulate!(secondOrder2, Tsit5(), stopTime=1.0, merge=Map(sys=Map(x=[1.0,1.0])),
          logParameters=true, logStates=true,
          requiredFinalStates = [1.0000295203337868, 0.0022367686372974493])
plot(secondOrder2, [("sys.u", "sys.y"); "sys.x"], figure = 3)
end
