module TestElectrical

using Modia
using Modia.Electric

# Desired:
#   using ModiaMath: plot
#   using Test
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Modia.Test
end

result = simulate(Resistor, 1)
@test result["i"][end] == 0.0
@test result["v"][end] == 0.0


@model ParallelResistors begin
    R1 = Resistor(R=1, p=Pin(v=Float(start=0.0)), n=Pin(v=Float(start=0.0)))
    R2 = Resistor(R=1)
    @equations begin
        connect(R1.p, R2.p)
        connect(R1.n, R2.n)
    end
end 

result = simulate(ParallelResistors, 1)
@test result["R1.i"][end] == 0.0
@test result["R1.v"][end] == 0.0

# -----------------------------------

# result = simulate(Capacitor, 1)

@model ParallelCapacitors begin
    C1 = Capacitor(C=2, v=Float(start=1.0))
    C2 = Capacitor(C=1, v=Float(state=false))
    ground = Ground()
    @equations begin
        connect(C1.p, C2.p)
        connect(C1.n, C2.n)
        connect(C1.n, ground.p)
    end
end 

result = simulate(ParallelCapacitors, 1)
@test result["C2.v"][end] == 1.0


# -----------------------------------

# result = simulate(Inductor, 1)

@model InductorsInSeries begin
    V = ConstantVoltage(V=10)
    L1 = Inductor(L=1, i=Current(start=0.0))
    R1 = Resistor(R=100) # , n=Pin(v=Float(start=0.0)))
    R2 = Resistor(R=100) # , n=Pin(v=Float(start=0.0)))
    L2 = Inductor(L=1, i=Current(start=0.0))
    @equations begin
        connect(V.p, L1.p)
        connect(L1.n, R1.p)
        connect(L1.n, R2.p)
        connect(R1.n, L2.p)
        connect(R2.n, L2.p)
        connect(L2.n, V.n)
    end
end 

# result = simulate(InductorsInSeries, 1)

# -----------------------------------

@model ParallelCapacitorCircuit begin
    R = Resistor(R=100.0)
    C = Capacitor(C=2.5E-2, v=Float(start=0.0))
    C2 = Capacitor(C=2.5E-2, v=Float(state=false))
    V = StepVoltage(startTime=10, V=100)
    #  Vt=ConstantVoltage(V=10.0)
    ground = Ground()
    @equations begin
        connect(V.n, ground.p)
        # connect(Vt.n, ground.p)
        connect(V.p, R.p)
        connect(R.n, C.p)
        connect(C.n, V.n)
        connect(C.n, C2.n)
        connect(C.p, C2.p)
    end
end 

# simulate(ParallelCapacitorCircuit, 1, useKinsol=true, removeSingularities=false, logTranslation=true)
# checkSimulation(ParallelCapacitorCircuit, 1, "C1.v", 1.0, useKinsol=true, removeSingularities=false)

end
