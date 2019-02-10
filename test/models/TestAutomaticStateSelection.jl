module TestAutomaticStateSelection

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

@model TwoConnectedInertias begin
    J1 = 2.0
    J2 = 3.0
    w1 = Float(size=(), start=1.0)
    w2 = Float(size=(), start=1.0)
    tau = Float(size=())
    t   = Float(size=())
    @equations begin
        der(t) = 1.0
        J1*der(w1) = sin(t) - tau
        J2*der(w2) = tau
        w1 = w2
    end
end 
result = simulate(TwoConnectedInertias, 3.0; logTranslation=true, logSimulation=true)
plot(result, ("w1", "w2"), figure=1)


@model TwoInertiasConnectedViaIdealGear begin
    J1 = 2.0
    J2 = 3.0
    ratio = 4.0
    w1 = Float(size=(), start=1.0)
    w2 = Float(size=(), start=1.0)
    tau = Float(size=())
    t   = Float(size=())
    @equations begin
        der(t) = 1.0
        J1*der(w1) = sin(t) - ratio*tau
        J2*der(w2) = tau
        w1 = ratio*w2
    end
end 
result = simulate(TwoInertiasConnectedViaIdealGear, 3.0; logTranslation=true, logSimulation=true)
plot(result, ("w1", "w2"), figure=2)



# Translation fails if tearing=true (successful with tearing=false)
# The reason is that "vj in Gsolvable[eq]" in Tearing.jl gives a wrong array access.
# It might be that the variables removed by RemoveSingularity are not correctly provided in Gsolvable.
#
@model ParallelCapacitors1 begin
    C1 = Capacitor(C=1.1, v=Float(start=1.0))
    C2 = Capacitor(C=2.2, v=Float(state=false))
    ground = Ground()
    @equations begin
        connect(C1.p, C2.p)
        connect(C1.n, C2.n)
        connect(C1.n, ground.p)
    end
end 

result = simulate(ParallelCapacitors1, 1; logTranslation=true, logSimulation=true, tearing=false, removeSingularities=true)
plot(result, ("C1.v", "C2.v"), figure=3)

 
@model ParallelCapacitors2 begin
    C1 = Capacitor(C=1.1, v=Float(start=1.0))
    C2 = Capacitor(C=2.2, v=Float(state=1.0))
    ground = Ground()
    @equations begin
        connect(C1.p, C2.p)
        connect(C1.n, C2.n)
        connect(C1.n, ground.p)
    end
end 

result = simulate(ParallelCapacitors2, 1; logTranslation=true, logSimulation=true, tearing=true, removeSingularities=true)
plot(result, ("C1.v", "C2.v"), figure=4)
end