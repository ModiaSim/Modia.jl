module Slides

using Modia

@define FirstOrder = Model(
    T = 0.5u"s",
    x = Var(init=0.3),
    equations = :[
        u = 2
        T * der(x) + x = u]
)

@define Pin = Model(v = potential, i = flow)

@define OnePort = Model(
    p = Pin,
    n = Pin,
    equations = :[
        0 = p.i + n.i
        v = p.v - n.v
        i = p.i ])

@define Resistor = OnePort(R = missing,
                           equations = :[ R*i = v ] )
                           
@define Capacitor = OnePort(C = missing,
                            v = Var(init=0.0u"V"), 
                            equations = :[ C*der(v) = i ] )
                            
@define ConstantVoltage = OnePort(V = 1.0u"V", equations = :[ v = V ] )

@define Circuit = Model(
    R = Resistor(R=0.5u"Ω"),
    C = Capacitor(C=2.0u"F"),
    V = ConstantVoltage(V=10.0u"V"),
    equations = :[
        connect(V.p, R.p)
        connect(R.n, C.p) 
        connect(C.n, V.n)])

end