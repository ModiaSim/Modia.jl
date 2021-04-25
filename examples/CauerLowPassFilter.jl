module CauerLowPassFilterModel

using Modia, ModiaPlot

l1 = 1.304
l2 = 0.8586
c1 = 1.072*1u"F"
c2 = 1/(1.704992^2*l1)*1u"F"
c3 = 1.682*1u"F"
c4 = 1/(1.179945^2*l2)*1u"F"
c5 = 0.7262*1u"F"

CauerLowPassOPVWithoutNodes(i=1) = Model(
  # Cauer low pass filter with operational amplifiers

  C1 = Modia.Capacitor | Map(C=c1 + c2), 
  C2 = Modia.Capacitor | Map(C=c2, v=Var(init=nothing)),
  C3 = Modia.Capacitor | Map(C=l1*1u"F"), 
  C4 = Modia.Capacitor | Map(C=c4*i, v=Var(init=nothing)),
  C5 = Modia.Capacitor | Map(C=c2, v=Var(init=nothing)), 
  R1 = Modia.Resistor | Map(R=1u"Ω"),
  R2 = Modia.Resistor | Map(R=1u"Ω"),
  R3 = Modia.Resistor | Map(R=1u"Ω"),
  Op1 = Modia.IdealOpAmp3Pin,
  G = Modia.Ground,
  R4 = Modia.Resistor | Map(R=-1u"Ω"),
  R5 = Modia.Resistor | Map(R=-1u"Ω"),
  Op2 = Modia.IdealOpAmp3Pin,
  Op3 = Modia.IdealOpAmp3Pin,
  G1 = Modia.Ground,
  R6 = Modia.Resistor | Map(R=1u"Ω"),
  R7 = Modia.Resistor | Map(R=1u"Ω"),
  C6 = Modia.Capacitor | Map(C=c2 + c3 + c4, v=Var(init=nothing)),
  R8 = Modia.Resistor | Map(R=-1u"Ω"),
  R9 = Modia.Resistor | Map(R=-1u"Ω"),
  R10 = Modia.Resistor | Map(R=1u"Ω"),
  Op4 = Modia.IdealOpAmp3Pin,
  Op5 = Modia.IdealOpAmp3Pin,
  C7 = Modia.Capacitor | Map(C=l2*1u"F"), 
  C8 = Modia.Capacitor | Map(C=c4, i=Var(start=0u"A")), 
  C9 = Modia.Capacitor | Map(C=c4 + c5), 
  R11 = Modia.Resistor | Map(R=1u"Ω"),

  G2 = Modia.Ground,
  G3 = Modia.Ground,
  G4 = Modia.Ground,
  V = Modia.ConstantVoltage | Map(V=1.0u"V"),
  Ground1 = Modia.Ground,

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
    filters = [CauerLowPassOPVWithoutNodes(0.1*i) for i in 1:10]
)

model = @instantiateModel(Filters, logDetails=false, logTiming=true, unitless=true)

println("Simulate")
@time simulate!(model, Tsit5(), stopTime = 60, requiredFinalStates = 
    [-0.5000732186007855, -0.5002239998029879, 0.4996923410661849, -0.4996706636198064, -0.5000929741546876, -0.5001255383344897, -0.500147742207576, 0.49981113006573824, -0.4996196061432069, -0.5001392063603706, -0.5001815616147427, -0.5000459266587362, 0.49996685095395915, -0.49958392015463465, -0.5001886259176451, -0.5002389821958331, -0.4999173001411904, 0.5001611780149184, -0.49957016523210945, -0.5002393712119582, -0.5002944833191254, -0.49976183273648583, 0.5003940588092074, -0.49958591550449066, -0.5002887040151874, -0.500343478728503, -0.4995812991530469, 0.5006629055320556, -0.49963969356758026, -0.5003327561376948, -0.50037980690654, -0.49938004306850625, 0.5009615326279456, -0.4997407922705924, -0.5003662177501006, -0.500395375828166, -0.4991659724894459, 0.5012787730322307, -0.4998989295998676, -0.5003819572251011, -0.5003797586125351, -0.49895184523346914, 0.5015966895590688, -0.5001236602409048, -0.5003705619890064, -0.5003197458974357, -0.49875691633150787, 0.5018882791340892, -0.5004234377294563, -0.5003197904255381])
plot(model, ["filters_1.C9.v", "filters_2.C9.v", "filters_3.C9.v"])

println("Total")
end
