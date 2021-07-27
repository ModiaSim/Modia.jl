module FilterCircuit

using Modia
@usingModiaPlot

setLogMerge(false)

include("$(Modia.modelsPath)/Electric.jl")

Filter = Model(
    R = Resistor | Map(R=0.5u"Ω"),
    C = Capacitor | Map(C=2.0u"F", v=Var(init=0.1u"V")),
    V = ConstantVoltage | Map(V=10.0u"V"),
    ground = Ground,
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n, ground.p)
    ]
)

model = @instantiateModel(Filter)
@time simulate!(model, Tsit5(), stopTime = 10, requiredFinalStates=[9.999550454584188]) 
plot(model, [("R.v", "C.v")])

println("Simulate once more with different R.R")
@time simulate!(model, Tsit5(), stopTime = 10, merge = Map(R = Map(R = 5u"Ω")), requiredFinalStates = [6.3579935215716095]) 
plot(model, [("R.v", "C.v")])

# setLogMerge(true)
println("Filter without ground and parameter propagation")
Filter2 = Model(
    r = 2.0u"Ω",
    c = 1.0u"F",
    v = 10u"V",
    R = Resistor | Map(R=:r),
    C = Capacitor | Map(C=:c),
    V = ConstantVoltage | Map(V=:v),
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n)
    ]
)

# @showModel(Filter2)

model = @instantiateModel(Filter2)
simulate!(model, Tsit5(), stopTime = 10, requiredFinalStates = [9.932620374719848]) 
plot(model, [("R.v", "C.v")])


println("Voltage divider by redeclaring capacitor to resistor")
Cpar = Map(C = 5.0u"F")

TwoFilters = Model( f1 = Filter | Map( R = Map(R = 10.0u"Ω"), C = Cpar), f2 = Filter) 

#setLogMerge(true)

VoltageDividerAndFilter = TwoFilters | Map(f1 = Map(C = redeclare | Resistor | (R = 2.0u"Ω", v = Map(start = 0u"V"))))

model = @instantiateModel(VoltageDividerAndFilter, unitless=true, logCalculations=false, evaluateParameters=true)
simulate!(model, Tsit5(), stopTime = 10) 
plot(model, [("f1.R.v", "f1.C.v"), ("f2.R.v", "f2.C.v")])

n=10
println("Build array of $n filters")
Filters = Model(
    filters = [Filter | Map( R = Map(R = (10.0+5*i)*u"Ω")) for i in 1:n]
)

setLogMerge(false)

model = @instantiateModel(Filters)

simulate!(model, Tsit5(), stopTime = 10, requiredFinalStates = 
[2.9063400246452358, 2.2898722474751163, 1.8945655444974656, 1.6198309235728083, 1.4179087924692246, 1.2632806644107324, 1.1410907635368692, 1.0421095614435398, 0.9603029088053439, 0.8915602951695468])
plot(model, [("filters_1.R.v", "filters_1.C.v"), ("filters_2.R.v", "filters_2.C.v")])

end