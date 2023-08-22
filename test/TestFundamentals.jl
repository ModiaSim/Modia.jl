module TestFundamentals

using Modia
@usingModiaPlot

setLogMerge(false)

include("../models/Electric.jl")

#=
@define Mod1 = Model(
    R1 = Resistor,
    R2 = Resistor | Map(R=10),
    R3 = Resistor(),
    R4 = Resistor(R=20),
    R5 = Resistor(R=10), # <: OnePort(v=0), # replaceable Resistor R5(R=10) constrainedby OnePort(v=0)
    CompType = Par(Resistor),
    R6 = CompType(),
    R7 = CompType(), # (v.init=1),
    equations = :[connect(R1.p, R2.n), connect(R2.p, R3.n), connect(R3.p, R4.n), 
    connect(R4.p, R5.n), connect(R5.p, R6.n), connect(R6.p, R7.n), connect(R7.p, R1.n)]
)

model = @instantiateModel(Mod1, logModel=true, aliasReduction=true, log=true, unitless=true, logCode=true)
simulate!(model, stopTime=1)
plot(model, ["R1.v", "R5.v"])

n = 10

@define Mod2 = Model(
#    p = Par(5),
#    q = Par(sin(p)),
#    u = Var(),
#    v = Var(cos(u)+n),
#    R = Resistor(R=p),   # R = Resistor | Map(R=10)
##    m1 = Mod1(R1(R=100)),   # m1 = Mod1(R1=Map(R=100))
#    m2 = Mod1(R1 = Capacitor),   # m2 = Mod1(R1=redeclare | Capacitor)
##    m2 = Mod1(R1 = Capacitor(C=5)),   # m2 = Mod1(R1=redeclare | Capacitor | Map(C=5))
#    MyResistor = Resistor | Model(equations=:[P=v*i]),   # extends Resistor
    m3 = Mod1A(CompType = Capacitor),   # Mod1 m3(redeclare model Comp=Resistor)

    m4 = [Resistor(R=i) for i in 1:5],
    m5 = if q > 0.5; Capacitor end

)

@define ModA = Model(
    CompType = Par(Resistor),
    comp = CompType()
)

@define ModB = Model(
    m = ModA(CompType = Capacitor),   # ModA m(redeclare model CompType=Capacitor)
)

@instantiateModel(ModB, logModel=true, log=true)

@define Circuit = Model(
    DeviceType = Par(Resistor),
    device = DeviceType(),
    current = ConstantCurrent(I=3),
    equations = :[
        connect(current.p, device.p)
        connect(current.n, device.n)
#        current.n.v = 0
    ]
)

@define CapacitorCircuit = Model(
    c = Circuit(DeviceType = Capacitor(C=2))   # ModA m(redeclare model CompType=Capacitor)
)

model = @instantiateModel(CapacitorCircuit, logModel=true, aliasReduction=true, log=true, unitless=true, logCode=true)
simulate!(model, stopTime=1)
plot(model, ["c.device.v"])

=#


@define RCCircuit = Model(
    r = Resistor(R=1),
    c = Capacitor(C=1),
    source = ConstantVoltage(V=1),
    equations = :[
        connect(source.p, r.p)
        connect(r.n, c.p)
        connect(source.n, c.n)
    ]
)

@define RCCircuits = Model(  # From slides
    rc1 = RCCircuit(r(R=4), c(C=10)),
    rc2 = RCCircuit(c(C=10, v(init=5)), source(V=-10))
)

@define RCCircuits2 = Model(
    rc1 = RCCircuit,
    rc2 = RCCircuit(),
    rc3 = RCCircuit | Map(r = Map(R = 2)) | Map(c = Map(C = 0.5)),
    rc4 = RCCircuit(r(R=2), c(C=0.5)), # Hierarchical modifiers (preferred)
#    rc5 = RCCircuit(r(R=2)) | Map(c = Map(C = 0.5)), # Mixed hierarchical and merge
    rc6 = RCCircuit(r.R=2, c.C=0.5),   # Dot notation
    rc7 = RCCircuit(r.R=2, c(C=0.5)),  # Mixed

    rc8 = RCCircuit(c(v(init=1)), source.V=-1), # Modifier for variable attribute
    rc9 = RCCircuit(c.v.init=1, source.V=-1), 
    rc10 = RCCircuit(c(v.init=1), source.V=-1), # Mixed
)

model = @instantiateModel(RCCircuits2, logModel=true, aliasReduction=true, log=true, unitless=true, logCode=true)
simulate!(model, stopTime=10)
plot(model, ["rc1.c.v", "rc2.c.v", "rc3.c.v", "rc4.c.v", "rc5.c.v", "rc6.c.v", "rc7.c.v", "rc8.c.v", "rc9.c.v", "rc10.c.v"])


@define RCCircuits3 = Model(
    rc3 = RCCircuit(r(R=2), c(C=0.5)), # Hierarchical modifiers
    rc3a = RCCircuit(r(R=2, ERR=5), r.R=10), # Multiple modifications, 
)

model = @instantiateModel(RCCircuits3, logModel=true, aliasReduction=true, log=true, unitless=true, logCode=true)
simulate!(model, stopTime=10)
plot(model, ["rc3.c.v", "rc3a.c.v"], figure=2)



@define Circuit = Model(
    DeviceType = Par(Resistor),
    device1 = DeviceType(),
    device2 = DeviceType(),
    source = ConstantVoltage(V=10),
    equations = :[
        connect(source.p, device1.p)
        connect(device1.n, device2.p)
        connect(source.n, device2.n)
    ]
)

#=
@define AlteredCircuit = Circuit(DeviceType = Resistor(R=2), device1(R=40))
@define AlteredCircuit2 = Circuit(device1(R=40))

@define CapacitorCircuit = Circuit(DeviceType = Capacitor, source=ConstantCurrent)

@define MixedCircuits1 = Model(
    r1 = Circuit(),  
    r2 = Circuit(DeviceType(R=2), device1(R=8)), 
    r3 = Circuit(DeviceType = Resistor(R=2), device1(R=8)),  

    c1 = Circuit(DeviceType = Capacitor(), source=ConstantCurrent(I=10)),
    c2 = Circuit(DeviceType = Capacitor(C=2), source=ConstantCurrent(I=10)),
    c3 = Circuit(DeviceType = Capacitor(C=2), device2(C=8), source=ConstantCurrent(I=10)),

    rc = Circuit(DeviceType = Capacitor(), device1=Resistor(R=2), device2(C=5), source(V=10))  

    #=
    c1 = CapacitorCircuit(source(I=10)),
    c2 = CapacitorCircuit(DeviceType = Capacitor(C=2), source(I=10)),
    c3 = CapacitorCircuit(DeviceType = Capacitor(C=3), device2(C=30), source(I=10))  
=#
    )
=#

@define MixedCircuits = Model(
    c = Circuit(DeviceType = Capacitor, device2(C=8), source=ConstantCurrent(I=10)),
    rc = Circuit(DeviceType = Capacitor, device1=Resistor(R=2), device2(C=5), source(V=10)), 
    rc1 = Circuit(DeviceType = Capacitor(C=5), device1=Resistor(R=2), source(V=10)), 
    rc2 = Circuit(device1=Capacitor(C=5), device2=Resistor(R=2), source(V=10))  
 
)

#=
@define GenericCircuits = Model(
    DevType = Capacitor,
    mc = MixedCircuits(rc(DeviceType=DevType))
)
=#

#=
@define CCCircuit = Model(
    c1 = Capacitor(C=2),
    c2 = Capacitor(C=4),
    source = ConstantCurrent(I=10),
    equations = :[
        connect(source.p, c1.p)
        connect(c1.n, c2.p)
        connect(source.n, c2.n)
    ]
)
=#

#=
model = @instantiateModel(MixedCircuits, logModel=true, aliasReduction=true, log=true, unitless=true, logCode=true)
simulate!(model, stopTime=100)
plot(model, ["c.device1.v", "c.device2.v"])
#plot(model, ["r3.device1.v", "r3.device2.v", "c3.device1.v", "c3.device2.v", "rc.device1.v", "rc.device2.v"])
#plot(model, ["rc.device1.v", "rc.device2.v"])
# plot(model, ["c1.v", "c2.v"])
=#

#=
@define Dummy = Model(
    DeviceType = Par(Resistor),
)

model = @instantiateModel(Dummy, logModel=true, aliasReduction=true, log=true, unitless=true, logCode=true)
simulate!(model, stopTime=100)
=#

end