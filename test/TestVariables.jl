module TestVariables

using Modia
@usingModiaPlot
using Modia.Measurements
using Modia.StaticArrays
using Modia.Test

@testset "... Test variable declarations" begin

v1 = input | output | flow | Var(min=0, nominal=10, init=5)
@test isequal(v1, Var(input = true, output = true, flow = true, min = 0, nominal = 10, init = 5))

v2 = parameter | 5 | Var(info="v2 info")
@test isequal(v2, Var(parameter = true, value = 5, info = "v2 info"))

v3 = "string value" | parameter | info"v3 info"
@test v3 == Var(parameter = true, value = "string value", info = "v3 info")

v4 = parameter | 5.5u"m/s" | Var(min=0)
@test isequal(v4, Var(parameter = true, value = 5.5u"m*s^-1", min = 0))

v5 = parameter | Var(5, min=0)
@test isequal(v5, Var(parameter = true, min = 0, value = 5))

Pin = Model( v = potential | Var(nominal=10), i = flow )
@test isequal(Pin, Model(v = Var(potential = true, nominal = 10), i = Var(flow = true)))

v6 = Var(2.0)
@test isequal(v6, Var(value = 2.0))

v7 = Var(2.0 ± 0.1)
@test isequal(v7, Var(value = 2.0 ± 0.1))

v8 = Var(SVector(1.0, 2.0, 3.0))
@test isequal(v8, Var(value = [1.0, 2.0, 3.0]))

v9 = v4 | Var(parameter=false)
@test isequal(v9, Var(parameter = false, value = 5.5u"m*s^-1", min = 0))

v10 = v4 | Var(parameter=nothing)
@test isequal(v10, Var(value = 5.5u"m*s^-1", min = 0))

v11 = v4 | interval(1,10)
@test isequal(v11, Var(parameter = true, value = 5.5u"m*s^-1", min = 1, max = 10))

TestVar1 = Model(
    p = parameter | -1 | info"A parameter",
    s = parameter | "string value" | info"A string parameter", # Not used
    x = Var(init=0),
    y = Var(start=0), # Not used
    equations = :[
        der(x) = p*x + 1
    ]
)

# @showModel TestVar1

model = @instantiateModel(TestVar1, logCode=true)
simulate!(model, stopTime= 1, requiredFinalStates = [1-exp(-1)]) 
plot(model, "x")

simulate!(model, stopTime= 1, merge=Map(p=-2, x=2), requiredFinalStates = [2-1.5*(1-exp(-2))]) 
plot(model, "x")

end

end