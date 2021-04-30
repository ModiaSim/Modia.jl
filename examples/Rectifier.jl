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
        Ron   = 1e-4,
        Goff  = 1e-4,
        s = Var(start = 0.0),
        equations = :[
            closed = positive(s)   # closed = s > 0
            v = s*(closed ? Ron : 1)
            i = s*(closed ? 1   : Goff)
        ]
    )


Rectifier = Model(
    R = Resistor | Map(R=1.0u"Ω"),
    C = Capacitor | Map(C=1.0u"F", v=Var(init=0.0u"V")),
    D = IdealDiode,
    V = SineVoltage | Map(V=5.0u"V", f=1.5u"Hz"),
    ground = Ground,
    connect = :[
      (V.p, D.p)
      (D.n, R.p)
      (C.p, R.p)
      (ground.p, R.n, C.n, V.n)
    ]
)

model = @instantiateModel(Rectifier, logCode=true)
@time simulate!(model, Tsit5(), stopTime = 3, nz = 1)
plot(model, ("R.v", "D.v", "V.v", "C.v"))

end