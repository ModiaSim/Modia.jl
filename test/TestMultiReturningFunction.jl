module TestMultiReturningFunction

using Modia
@usingModiaPlot

include("$(Modia.path)/models/Blocks.jl")
include("$(Modia.path)/models/Electric.jl")
include("$(Modia.path)/models/Rotational.jl")


function twoCoupledInertias(J1,J2,derw1,derw2,tau0,tau1,tau2)::Array{Float64,1}
    r1 = J1*derw1 - tau0 + tau1
    r2 = J2*derw2 - tau1 + tau2
    return [r1,r2]
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
    
    # Same model with components from the Rotational.jl library
    inertia1 = Inertia | Map(J=:J1),
    inertia2 = Inertia | Map(J=:J2),
    inertia3 = Inertia | Map(J=:J3),
    connect = :[        
        (inertia1.flange_b, inertia2.flange_a)
        (inertia2.flange_b, inertia3.flange_a)
    ],
       
    equations = :[
        tau0 = sin(time/u"s")    
        w1 = der(phi1)
        w2 = der(phi2)
        w3 = der(phi3)
        phi2 = phi1
        phi3 = phi2
        (0,0) = twoCoupledInertias(J1,J2,der(w1),der(w2),tau0,tau1,tau2)
        J3*der(w3) = tau2
        inertia1.flange_a.tau = tau0
        inertia3.flange_b.tau = 0.0
    ]
)
    
threeCoupledInertias = @instantiateModel(ThreeCoupledInertias, unitless=true, log=false, logDetails=false, logCode=true, logStateSelection=false)

simulate!(threeCoupledInertias, stopTime = 2.0, log=true, 
          requiredFinalStates=[0.3029735305821086, 0.3933746028781303, 0.39337460287813036, 0.3029735305821086])

plot(threeCoupledInertias, [("phi1", "phi2", "phi3", "inertia1.phi"), 
                            ("w1","w2","w3", "inertia1.w"), 
                            ("tau0"),
                            ("tau1", "inertia2.flange_a.tau"),                            
                            ("tau2", "inertia3.flange_a.tau")])

end
