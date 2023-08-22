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


include("HeatTransfer/InsulatedRod2.jl")

"""
    insulatedRod = InsulatedRod2(; L, A, rho=7500.0u"kg/m^3", lambda=74.0u"W/(m*K)", c=450.0u"J/(kg*K), T0=293.15u"K", nT=1") 
    
Generate a Model(..) instance of an insulated rod with length `L` and cross sectional area `A` that
models 1D heat transfer from port_a to port_b. The rod is discretized with `nT` internal temperature nodes that are initialized with `T0`.
This is an acausal built-in component where the code size does not depend on `nT` and
`nT` can still be changed after model code is generated and compiled.

For more details see Appendix B1 of [DOI: 10.3390/electronics12030500](https://doi.org/10.3390/electronics12030500).

# Optional Arguments:
- `rho`: density of rod material
- `lambda`: thermal conductivity of rod material
- `c`: specific heat capacity of rod material
- `T0`: initial temperature of internal temperature nodes
- `nT`: number of temperature nodes (`nT > 0`).

# Internal variables that can be inquired/plotted:
- `T`: Vector of temperatures at the internal temperature nodes
- `der(T)': Vector of derivatives of `T`.
"""
InsulatedRod2 = Model(; 
    _buildFunction       = Par(functionName = :(build_InsulatedRod2!)),        # Called once in @instantiateModel(..) before getDerivatives!(..) is generated 
    _initSegmentFunction = Par(functionName = :(initSegment_InsulatedRod2!)),  # Called once before initialization of a new simulation segment
    L      = 1.0u"m",                # Length of rod
    A      = 0.0004u"m^2",           # Rod area
    rho    = 7500.0u"kg/m^3",        # Density of rod material
    lambda = 74.0u"W/(m*K)",         # Thermal conductivity of rod material
    c      = 450.0u"J/(kg*K)",       # Specific heat capacity of rod material
    T0     = 293.15u"K",             # Initial temperature of internal temperature nodes
    nT     = 1,                      # Number of temperature nodes (nT > 0).
    port_a  = HeatPort,              # heat port on left side
    port_b  = HeatPort               # heat port on right side
)
