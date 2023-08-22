module PendulumSimulation

using Modia
using Modia.Measurements
@usingModiaPlot


Pendulum = Model(
   L = 0.8u"m",
   m = 1.0u"kg",
   d = 0.5u"N*m*s/rad",
   g = 9.81u"m/s^2",
   phi = Var(init = 1.57*u"rad"),
   w   = Var(init = 0u"rad/s"),
   equations = :[
          w = der(phi)
        0.0 = m*L^2*der(w) + d*w + m*g*L*sin(phi)
          r = [L*cos(phi), -L*sin(phi)]
   ]
)


pendulum = @instantiateModel(Pendulum)
simulate!(pendulum, Tsit5(), stopTime = 10.0u"s", log=true)
plot(pendulum, [("phi", "w"); "r"])

PendulumWithUncertainties = Pendulum | Map(L = (0.8 ± 0.2)u"m",
                                           m = (1.0 ± 0.2)u"kg",
                                           d = (0.5 ± 0.2)u"N*m*s/rad")

pendulum2 =  @instantiateModel(PendulumWithUncertainties,
                               FloatType = Measurement{Float64})

simulate!(pendulum2, Tsit5(), stopTime = 10.0u"s")
plot(pendulum2, [("phi", "w"); "r"], figure = 2)


# Linearize
println("\n... Numerically linearize at stopTime = 10 with Float64 and Double64:")
(A_10, x_10) = linearize!(pendulum2, stopTime=10) 

#= DoubleFloats is not necessarily defined in the environment
using Modia.DoubleFloats
pendulum3 = InstantiatedModel{Measurement{Double64}}(pendulum2)
(A_10_Double64, x_10_Double64) = linearize!(pendulum3, stopTime=10) 
=#

xNames = get_xNames(pendulum2)
@show xNames
println(IOContext(stdout, :error_digits=>15), "A_10 = ", A_10, ", x_10 = ", x_10)
#println(IOContext(stdout, :error_digits=>15), "A_10_Double64 = ", A_10_Double64, ", x_10_Double64 = ", x_10_Double64)

end