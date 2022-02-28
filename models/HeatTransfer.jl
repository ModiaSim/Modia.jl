#=
Modia library with 1D heat transfer component models (inspired from Modelica Standard Library).

Developer: Martin Otter, DLR-SR  
Copyright 2021, DLR Institute of System Dynamics and Control
License: MIT (expat)
=#

using Modia

HeatPort = Model( T = potential,   # Absolute temperature
                  Q_flow = flow )  # Heat flow into the component

FixedTemperature = Model(
    T    = 293.15u"K",
    port = HeatPort,
    equations = :[port.T = T]
)


FixedHeatFlow = Model(
    Q_flow = 0.0u"W",   # Fixed heat flow into the connected component
    port   = HeatPort,
    equations = :[port.Q_flow = -Q_flow]
)


HeatCapacitor = Model(
    C = 1.0u"J/K",
    port = HeatPort,
    T = Var(init = 293.15u"K"),
    equations = :[T = port.T,
                  der(T) = port.Q_flow/C]
)


ThermalConductor = Model(
    G = 1.0u"W/K",
    port_a = HeatPort,
    port_b = HeatPort,
    equations = :[dT = port_a.T - port_b.T
                  0  = port_a.Q_flow + port_b.Q_flow
                  port_a.Q_flow = G*dT]
)



# Fully insulated rod with 1D heat transfer and port_a/port_b on left/right side

T_grad1(T,Ta,dx,i) = i == 1         ? (Ta - T[1])/(dx/2) : (T[i-1] - T[i]  )/dx
T_grad2(T,Tb,dx,i) = i == length(T) ? (T[i] - Tb)/(dx/2) : (T[i]   - T[i+1])/dx

InsulatedRod = Model( 
    L       = 1.0u"m",            # Length of rod
    A       = 0.0004u"m^2",       # Rod area
    rho     = 7500.0u"kg/m^3",    # Density of rod material
    lambda  = 74.0u"W/(m*K)",     # Thermal conductivity of rod material
    c       = 450.0u"J/(kg*K)",   # Specific heat capacity of rod material
    port_a  = HeatPort,           # Heat port on left side
    port_b  = HeatPort,           # Heat port on right side
    T       = Var(init = fill(293.15u"K", 1)),  # Temperatures at the internal nodes
    equations = :[
        n   = length(T)
        dx  = L/n
        Ce  = c*rho*A*dx
        k   = lambda*A/Ce
        der(T) = k*[T_grad1(T,port_a.T,dx,i) - T_grad2(T,port_b.T,dx,i) for i in eachindex(T)]
        port_a.Q_flow = lambda*A*(port_a.T - T[1])/(dx/2)
        port_b.Q_flow = lambda*A*(port_b.T - T[n])/(dx/2)      
    ]
)

    
