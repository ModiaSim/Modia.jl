module TestModifierMode

using Modia
@usingModiaPlot

include("../models/Electric.jl")

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

model = @instantiateModel(RCCircuits, unitless=true)
simulate!(model, stopTime=10, requiredFinalStates=[0.22119922721253596, -4.48181591672524])
plot(model, ["rc1.c.v", "rc2.c.v"], figure=1)


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

model = @instantiateModel(RCCircuits2, unitless=true)
simulate!(model, stopTime=10, requiredFinalStates=[0.9999543756364957, 0.9999543756364957, 0.9999543756364957, 0.9999543756364957, 0.9999543756364957, 0.9999543756364957, -0.9999087910440049, -0.9999087910440049, -0.9999087910440049])
plot(model, ["rc1.c.v", "rc2.c.v", "rc3.c.v", "rc4.c.v",   "rc6.c.v", "rc7.c.v", "rc8.c.v", "rc9.c.v", "rc10.c.v"], figure=2)


@define RCCircuits3 = Model(
    rc3 = RCCircuit(r(R=2), c(C=0.5)), # Hierarchical modifiers
    rc3a = RCCircuit(r(R=2, ERR=5), r.R=10), # Multiple modifications, 
)

model = @instantiateModel(RCCircuits3, unitless=true)
simulate!(model, stopTime=10, requiredFinalStates=[0.9999565762239844, 0.632120613764388])
plot(model, ["rc3.c.v", "rc3a.c.v"], figure=3)


@define Circuit = Model(
    DeviceType = Resistor, # Generic model parameter
    device1 = DeviceType(),
    device2 = DeviceType(),
    source = ConstantVoltage(V=10),
    equations = :[
        connect(source.p, device1.p)
        connect(device1.n, device2.p)
        connect(source.n, device2.n)
    ]
)

@define MixedCircuits = Model(
    c1 = Circuit(),
    c2 = Circuit(device1(R=1), device2(R=4)),
    c3 = Circuit(DeviceType=Resistor(R=2), device2(R=8)),

#    cc1 = Circuit(device1=Capacitor(C=5), device2=Capacitor(C=20), source=ConstantCurrent(I=10)),  
    cc2 = Circuit(DeviceType = Capacitor(C=0.5), device2(C=2), source=ConstantCurrent(I=10)),

#    rc1 = Circuit(device1=Capacitor(C=5), device2=Resistor(R=2), source(V=10)),  
    rc2 = Circuit(DeviceType = Capacitor(C=5), device1=Resistor(R=2), source(V=10)), 
    rc3 = Circuit(DeviceType = Capacitor(), device1=Resistor(R=2), device2(C=5), source(V=10)), 
)

model = @instantiateModel(MixedCircuits, unitless=true)
simulate!(model, stopTime=100) #, requiredFinalStates=[-999.9999999999999, -124.99999999999999])
plot(model, ["c1.device1.v", "c1.device2.v"], figure=4)
plot(model, ["c2.device1.v", "c2.device2.v"], figure=5)
plot(model, ["c3.device1.v", "c3.device2.v"], figure=6)
plot(model, ["cc2.device1.v", "cc2.device2.v"], figure=7)
plot(model, ["rc2.device1.v", "rc2.device2.v"], figure=8)
plot(model, ["rc3.device1.v", "rc3.device2.v"], figure=9)

end