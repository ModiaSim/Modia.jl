module TestVariables

using TinyModia
using ModiaPlot
using Measurements
using StaticArrays
using Test

@testset "... Test variable declarations" begin

v1 = input | output | flow | Var(min=0, nominal=10, init=5)
@test v1 == (class = :Var, input = true, output = true, flow = true, min = 0, nominal = 10, init = 5)

v2 = parameter | 5 | Var(info="v2 info")
@test v2 == (class = :Var, parameter = true, value = 5, info = "v2 info")

v3 = "string value" | parameter | info"v3 info"
@test v3 == (class = :Var, parameter = true, value = "string value", info = "v3 info")

v4 = parameter | 5.5u"m/s" | Var(min=0)
@test v4 == (class = :Var, parameter = true, value = 5.5u"m*s^-1", min = 0)

v5 = parameter | Var(5, min=0)
@test v5 == (class = :Var, parameter = true, value = 5, min = 0)

Pin = Model( v = potential | Var(nominal=10), i = flow )
@test Pin == (class = :Model, v = (class = :Var, potential = true, nominal = 10), i = (class = :Var, flow = true,))

v6 = Var(2.0)
@test v6 == (class = :Var, value = 2.0,)

v7 = Var(2.0 ± 0.1)
@test v7 == (class = :Var, value = 2.0 ± 0.1,)

v8 = Var(SVector(1.0, 2.0, 3.0))
@test v8 == (class = :Var, value = [1.0, 2.0, 3.0],)

v9 = v4 | Var(parameter=false)
@test v9 == (class = :Var, parameter = false, value = 5.5u"m*s^-1", min = 0)

v10 = v4 | Var(parameter=nothing)
@test v10 == (class = :Var, value = 5.5u"m*s^-1", min = 0)

v11 = v4 | interval(1,10)
@test v11 == (class = :Var, parameter = true, value = 5.5u"m*s^-1", min = 1, max = 10)

TestVar1 = Model(
    p = parameter | -1 | info"A parameter",
    s = parameter | "string value" | info"A string parameter",
    x = Var(init=0.1),
    y = Var(start=0.1), # Not used
    equations = :[
        der(x) = p*x+1
    ]
)
#@showModel TestVar1
model = @instantiateModel(TestVar1, log=false, logCode=false)
simulate!(model, merge=Map(p=-2, x=0.2), requiredFinalStates = [0.4593994150057028]) # x.init is not changed
plot(model, "x")

end

end