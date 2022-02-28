module TestMultiReturningFunction4A

using Modia
using Modia.StaticArrays
@usingModiaPlot

mutable struct Mbs
    derw1
    derw2
end

function twoCoupledInertias(mbs::Mbs,J1,J2,tau0)
    tau1 = tau0 - J1*mbs.derw1
    tau2 = tau1 - J2*mbs.derw2
    return (tau1,tau2)
end


ThreeCoupledInertias = Model(
    J1 = 1.1,
    J2 = 1.2,
    J3 = 1.3,
    phi1 = Var(init=0.0),
    phi2 = Var(init=0.0),
    phi3 = Var(init=0.0),
    w1   = Var(init=0.0),
    w2   = Var(init=0.0),
    w3   = Var(init=0.0),
    tau1 = Var(start=0.0),
    tau2 = Var(start=0.0),
    equations = :[
        tau0 = sin(time/u"s")    
        w1 = der(phi1)
        w2 = der(phi2)
        w3 = der(phi3)
        phi2 = phi1
        phi3 = phi2
        mbs = Mbs(der(w1), der(w2))
        (tau1,tau2) = twoCoupledInertias(mbs,J1,J2,tau0)
        J3*der(w3) = tau2
    ]
)
    
threeCoupledInertias = @instantiateModel(ThreeCoupledInertias, unitless=true, log=false, logDetails=false, logCode=true, logStateSelection=false)
simulate!(threeCoupledInertias, stopTime = 2.0, log=true, requiredFinalStates=[0.3933746028781301, 0.3029735305821084])

plot(threeCoupledInertias, [("phi1", "phi2", "phi3"), 
                            ("w1","w2","w3"), 
                            ("tau0", "tau1", "tau2")])

end
