module CauerLowPassFilterModel

using TinyModia

using DifferentialEquations
using ModiaPlot

include("../models/Electric.jl")

l1 = 1.304
l2 = 0.8586
c1 = 1.072*1u"F"
c2 = 1/(1.704992^2*l1)*1u"F"
c3 = 1.682*1u"F"
c4 = 1/(1.179945^2*l2)*1u"F"
c5 = 0.7262*1u"F"

CauerLowPassOPVWithoutNodes(i=1) = Model(
  # Cauer low pass filter with operational amplifiers

  C1 = Capacitor | Map(C=c1 + c2), 
  C2 = Capacitor | Map(C=c2, v=Var(init=nothing)),
  C3 = Capacitor | Map(C=l1*1u"F"), 
  C4 = Capacitor | Map(C=c4*i, v=Var(init=nothing)),
  C5 = Capacitor | Map(C=c2, v=Var(init=nothing)), 
  R1 = Resistor | Map(R=1u"Ω"),
  R2 = Resistor | Map(R=1u"Ω"),
  R3 = Resistor | Map(R=1u"Ω"),
  Op1 = IdealOpAmp3Pin,
  G = Ground,
  R4 = Resistor | Map(R=-1u"Ω"),
  R5 = Resistor | Map(R=-1u"Ω"),
  Op2 = IdealOpAmp3Pin,
  Op3 = IdealOpAmp3Pin,
  G1 = Ground,
  R6 = Resistor | Map(R=1u"Ω"),
  R7 = Resistor | Map(R=1u"Ω"),
  C6 = Capacitor | Map(C=c2 + c3 + c4, v=Var(init=nothing)),
  R8 = Resistor | Map(R=-1u"Ω"),
  R9 = Resistor | Map(R=-1u"Ω"),
  R10 = Resistor | Map(R=1u"Ω"),
  Op4 = IdealOpAmp3Pin,
  Op5 = IdealOpAmp3Pin,
  C7 = Capacitor | Map(C=l2*1u"F"), 
  C8 = Capacitor | Map(C=c4, i=Var(start=0u"A")), 
  C9 = Capacitor | Map(C=c4 + c5), 
  R11 = Resistor | Map(R=1u"Ω"),

  G2 = Ground,
  G3 = Ground,
  G4 = Ground,
  V = ConstantVoltage | Map(V=1.0u"V"),
  Ground1 = Ground,

connect = :[ 
  (Op1.in_p, G.p) 
  (G1.p, Op2.in_p) 
  (R1.n, Op1.in_n, C2.n, R2.n, C1.p, R3.p) 
  (R3.n, C1.n, Op1.out, R4.p, C5.p) 
  (R4.n, Op2.in_n, C3.p, R5.n) 
  (C2.p, R5.p, Op3.out, C6.n, R9.p, C8.p) 
  (C3.n, Op2.out, R2.p, R7.p)
  (R7.n, Op3.in_n, C5.n, R6.n, C6.p, C4.n) 
  (C4.p, R8.p, R11.n, C9.n, Op5.out)
  (R9.n, R8.n, Op4.in_n, C7.p) 
  (R6.p, C7.n, Op4.out, R10.p)
  (G2.p, Op3.in_p) 
  (R11.p, C9.p, R10.n, Op5.in_n, C8.n) 
  (Op4.in_p, G3.p) 
  (Op5.in_p, G4.p) 
  (V.p, Ground1.p)            
  (V.n, R1.p)  ]
)

println("Build array of Cauer low pass filters")
@time Filters = Model(
    filters = [CauerLowPassOPVWithoutNodes(1) for i in 1:6]
)

model = @instantiateModel(Filters, logDetails=false, logTiming=false, unitless=true)

println("Simulate")
@time simulate!(model, Tsit5(), stopTime = 60)
plot(model, ("filters_1.C9.v", 
             "filters_2.C9.v",
             "filters_3.C9.v",
             "filters_4.C9.v",  
             "filters_5.C9.v",    
             "filters_6.C9.v",     
             ) )
end
