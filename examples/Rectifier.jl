module RectifierSimulation

using TinyModia
using DifferentialEquations
using ModiaPlot

setLogMerge(false)

include("../models/Electric.jl")


# Sinusoidal voltage source
SineVoltage = OnePort | Model( V = 1.0u"V", f = 1.0u"Hz", equations = :[ v = V*sin(2*3.14*f*time) ] )

# Ideal diode
IdealDiode = OnePort | Model(
        Ron   = 1e-5u"Ω",
        Goff  = 1e-5u"1/Ω",
        s = Var(start = 0.0),
        equations = :[              
            closed = positive(s)   # closed = s > 0  
            v = s*u"V"*(if closed; Ron*u"1/Ω" else 1         end)
            i = s*u"A"*(if closed; 1          else Goff*u"Ω" end)
        ]
    )


Rectifier1 = Model(
    R1 = Resistor | Map(R=1.0u"Ω"),
    R2 = Resistor | Map(R=100.0u"Ω"),    
    C = Capacitor | Map(C=0.1u"F", v=Var(init=0.0u"V")),
    D = IdealDiode,
    V = SineVoltage | Map(V=5.0u"V", f=1.5u"Hz"),
    ground = Ground,
    connect = :[
      (V.p , R1.p)
      (R1.n, D.p)
      (D.n , R2.p, C.p)
      (ground.p, R2.n, C.n, V.n)
    ]
)

rectifier1 = @instantiateModel(Rectifier1)
@time simulate!(rectifier1, Tsit5(), stopTime = 3, requiredFinalStates=[4.5665505086780565])
plot(rectifier1, [("V.v", "D.v", "C.v"), "D.i"], figure=1)


Rectifier2 = Model(
    V  = SineVoltage | Map(V=220.0u"V", f=50.0u"Hz"),
    R1 = Resistor    | Map(R=20.0u"Ω"),
    R2 = Resistor    | Map(R=500.0u"Ω"),    
    C  = Capacitor   | Map(C=1e-4u"F", v=Var(init=0.0u"V")),
    D1 = IdealDiode,
    D2 = IdealDiode,
    D3 = IdealDiode,
    D4 = IdealDiode,  
    ground = Ground,
    connect = :[
      (V.p , R1.p)
      (R1.n, D1.p, D3.n)
      (V.n , D2.p, D4.n, ground.p)
      (D1.n, D2.n, R2.p, C.p)
      (D3.p, D4.p, R2.n, C.n)
    ]
)

rectifier2 = @instantiateModel(Rectifier2, logExecution=false, logCode=false, unitless=true)
@time simulate!(rectifier2, Tsit5(), stopTime = 0.1, requiredFinalStates=[183.9899542497182])
plot(rectifier2, [("V.v", "C.v"), "V.i"], figure=2)


end