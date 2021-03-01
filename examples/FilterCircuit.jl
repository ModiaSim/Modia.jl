module FilterCircuit

using TinyModia

using DifferentialEquations
using ModiaPlot

setLogMerge(false)

include("../models/Electric.jl")

# res = Resistor()
# @showModel(res)

Filter = Model(
    R = Resistor | Map(R=0.5u"Ω"),
    C = Capacitor | Map(C=2.0u"F", init=Map(v=0.1u"V")),
    V = ConstantVoltage | Map(V=10.0u"V"),
    ground = Ground,
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n, ground.p)
    ]
)

Filter2 = Model(
    r = 1.0u"Ω",
    c = 1.0u"F",
    v = 1.0u"V",
    R = Resistor | Map(R=:(up.r)),
    C = Capacitor | Map(C=:(up.c)),
    V = ConstantVoltage | Map(V=:(up.v)),
    ground = Ground,
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n, ground.p)
    ]
)

Cpar = Map(C = 5.0u"F")

TwoFilters = Model( f1 = Filter | Map( R = Map(R = 10.0u"Ω"), C = Cpar), f2 = Filter) 

VoltageDividerAndFilter = TwoFilters | Map(f1 = Map(C = Redeclare | Resistor | (R = 20.0u"Ω", start = Map(v = 0u"V"))))

println("Build")
@time Filters = Model(
    filters = [Filter for i in 1:10]
)

setLogMerge(false)

model = @instantiateModel(Filters, logDetails=false, aliasReduction=false)

println("Simulate")
@time simulate!(model, Tsit5(), stopTime = 50, requiredFinalStates = 
    [9.999999294887072, 9.999999294887072, 9.999999294887072, 9.999999294887072, 9.999999294887072, 9.999999294887072, 9.999999294887072, 9.999999294887072, 9.999999294887072, 9.999999294887072])
plot(model, [("filters_1.R.v", "filters_1.C.v"), ("filters_2.R.v", "filters_2.C.v")])

end