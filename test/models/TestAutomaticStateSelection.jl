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

const logSimulation=true


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
result = simulate(TwoConnectedInertias, 3.0; logTranslation=true, logSimulation=logSimulation, tearing=true, removeSingularities=false, automaticStateSelection=true)
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
result = simulate(TwoInertiasConnectedViaIdealGear, 3.0; logTranslation=true, logSimulation=logSimulation, tearing=true, removeSingularities=true, automaticStateSelection=true)
plot(result, ("w1", "w2"), figure=2)

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

result = simulate(ParallelCapacitors1, 1; logTranslation=true, logSimulation=logSimulation, tearing=true, removeSingularities=true, automaticStateSelection=false)
plot(result, ("C1.v", "C2.v"), figure=3)
 
@model ParallelCapacitors2 begin
    C1 = Capacitor(C=1.1, v=Float(start=1.0))
    C2 = Capacitor(C=2.2, v=Float(start=1.0))
    ground = Ground()
    @equations begin
        connect(C1.p, C2.p)
        connect(C1.n, C2.n)
        connect(C1.n, ground.p)
    end
end 

result = simulate(ParallelCapacitors2, 1; logTranslation=true, logSimulation=logSimulation, tearing=true, removeSingularities=false, automaticStateSelection=true)
plot(result, ("C1.v", "C2.v"), figure=4)


@model ParallelCapacitors2b begin
    C1 = 1e-3
    C2 = 2e-3
    u1 = Float(size=(), start=1.0)
    u2 = Float(size=(), state=false)
    i1 = Float(size=())        
    v1 = Float(size=())        
    v0 = Float(size=())        
@equations begin
    C1*der(u1) = i1
    C2*der(u2) = -i1
    u1 = v1 - v0
    u2 = v1 - v0
    v0 = 0
    end 
end 
result = simulate(ParallelCapacitors2b , 1.0; logTranslation=true, logSimulation=logSimulation, tearing=true, removeSingularities=true, automaticStateSelection=false)
plot(result, ("u1", "u2", "i1", "v1"), figure=6)


# Does not compile
@model TwoInertiasConnectedViaIdealGearWithPositionConstraints begin
    J1 = 2.0
    J2 = 3.0
    ratio = 4.0
    phi1 = Float(size=(), start=1.0, fixed=true)
    phi2 = Float(size=(), start=1.0)
    w1 = Float(size=(), start=1.0, fixed=true)
    w2 = Float(size=(), start=1.0, fixed=true)
    tau = Float(size=())
    t   = Float(size=())
    @equations begin
        der(t) = 1.0
        w1 = der(phi1)
        w2 = der(phi2)
        J1*der(w1) = sin(t) - ratio*tau
        J2*der(w2) = tau
        phi1 = ratio*phi2
    end
end 
result = simulate(TwoInertiasConnectedViaIdealGearWithPositionConstraints, 3.0; logTranslation=true, logSimulation=logSimulation, tearing=true, removeSingularities=false, automaticStateSelection=true)
plot(result, [("phi1", "phi2"), ("w1", "w2")], figure=7)
println("after plot")


end