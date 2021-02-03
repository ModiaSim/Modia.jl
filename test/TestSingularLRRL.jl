module TestSingularLRRL

using TinyModia
using DifferentialEquations
using ModiaPlot


SingularLRRL = Model(
    v0   = 10,
    L1_L = 0.1,
    L2_L = 0.2,
    R1_R = 10,
    R2_R = 20,

    init = Map(R1_p_i = 0.1),
    equations =:[
        # Inductor 1
        L1_v = L1_p_v - L1_n_v
        0 = L1_p_i + L1_n_i
        L1_L*der(L1_p_i) = L1_v

        # Inductor 2
        L2_v = L2_p_v - L2_n_v
        0 = L2_p_i + L2_n_i
        L2_L*der(L2_p_i) = L2_v
    
        # Resistor 1
        R1_v = R1_p_v - R1_n_v
        0 = R1_p_i + R1_n_i
        R1_v = R1_R*R1_p_i
    
        # Resistor 2
        R2_v = R2_p_v - R2_n_v
        0 = R2_p_i + R2_n_i
        R2_v = R2_R*R2_p_i
        
        # Voltage source
        V_v = V_p_v - V_n_v
        0   = V_p_i + V_n_i
        
        # Connect equations
        #   connect(V.p, L1.p)
        #   connect(L1.n, R1.p)
        #   connect(L1.n, R2.p)
        #   connect(R1.n, L2.p)
        #   connect(R2.n, L2.p)
        #   connect(L2.n, V.n)
        V_p_v = L1_p_v
        0 = V_p_i + L1_p_i
        L1_n_v = R1_p_v
        L1_n_v = R2_p_v                
        0 = L1_n_i + R1_p_i + R2_p_i
        R1_n_v = L2_p_v
        R2_n_v = L2_p_v                
        0 = R1_n_i + R2_n_i + L2_p_i
        L2_n_v = V_n_v
        0 = L2_n_i + V_n_i #+ ground_v_i
        #V_n_v = 0
        V_v = v0 ]
)

singularLRRL = @instantiateModel(SingularLRRL)

simulate!(singularLRRL, Tsit5(), stopTime = 1.0)

plot(singularLRRL, [("L1_p_i", "L2_p_i", "R1_p_i", "R2_p_i"), ("L1_v", "L2_v", "R1_v", "R2_v")])

include("../models/Electric.jl")  

setLogMerge(false)

SingularLRRL2 = Model(
    R1 = Resistor | Map(R=10),
    R2 = Resistor | Map(R=20, init = Map(v=2)),
    L1 = Inductor | Map(L=0.1, init = Map(i=nothing)),
    L2 = Inductor | Map(L=0.2, init = Map(i=nothing)),
    V = ConstantVoltage | Map(V=10),
    connections = :[
        (V.p, L1.p)
        (L1.n, R1.p, R2.p)
        (L2.p, R1.n, R2.n)
        (L2.n, V.n) ]
)

singularLRRL2 = @instantiateModel(SingularLRRL2)

simulate!(singularLRRL2, Tsit5(), stopTime = 1.0)

plot(singularLRRL2, [("L1.p.i", "L2.p.i", "R1.p.i", "R2.p.i"), ("L1.v", "L2.v", "R1.v", "R2.v")])

SingularLRRL3 = Model(
    R1 = Resistor | Map(R=10),
    R2 = Resistor | Map(R=20), #| Map(init = Map(v=2)),
    R3 = Resistor | Map(R=30),
    R4 = Resistor | Map(R=40),
    L1 = Inductor | Map(L=0.1), # | Map(init = Map(i=nothing)),
    L2 = Inductor | Map(L=0.2), # | Map(init = Map(i=nothing)),
    V = ConstantVoltage | Map(V=10),
    ground = Ground,
    connections = :[
        (V.p, L1.p)
        (L1.n, R1.p, R2.p)
        (R1.n, R2.n, R3.p, R4.p)
        (L2.p, R3.n, R4.n)
        (L2.n, V.n, ground.p) ]
)

#=
singularLRRL3 = @instantiateModel(SingularLRRL3, unitless=true, log=true, logDetails=true)

simulate!(singularLRRL3, Tsit5(), stopTime = 1.0)

plot(singularLRRL3, [("L1.p.i", "L2.p.i", "R1.p.i", "R2.p.i"), ("L1.v", "L2.v", "R1.v", "R2.v")])
=#

end
