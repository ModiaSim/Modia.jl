module TestFilterCircuit2

using Modia
@usingModiaPlot

setLogMerge(false)

include("../models/Electric.jl")

FilterCircuit2 = Model(
    R  = Resistor  | Map(R=100u"Ω"),
    Ri = Resistor  | Map(R=10u"Ω"),
    C  = Capacitor | Map(C=2.5e-3u"F", v=Var(init=0.0u"V")),
    V  = ConstantVoltage | Map(V=10.0u"V"),
    ground = Ground,     
    connect = :[
      (V.p, Ri.n)
      (Ri.p, R.p)
      (R.n, C.p)
      (C.n, V.n, ground.p)
    ]
)

filterCircuit2 = @instantiateModel(FilterCircuit2, unitless=true, logCode=true)

simulate!(filterCircuit2, Tsit5(), stopTime = 1.0, requiredFinalStates =[9.736520087687524]) 
plot(filterCircuit2, ("V.v", "Ri.v", "R.v", "C.v"), figure=1)

end