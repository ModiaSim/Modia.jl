module TestFunctionWithSeveralReturnArguments

using ModiaLang
@usingModiaPlot


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
    equations = :[
        tau0 = sin(time/u"s")    
        w1 = der(phi1)
        w2 = der(phi2)
        w3 = der(phi3)
        phi2 = phi1
        phi3 = phi2
        (0,0) = twoCoupledInertias(J1,J2,der(w1),der(w2),tau0,tau1,tau2)
        J3*der(w3) = tau2
    ]
)
    
threeCoupledInertias = @instantiateModel(ThreeCoupledInertias, log=true, logDetails=true, logCode=true, logStateSelection=true)


simulate!(threeCoupledInertias, stopTime = 2.0)

plot(threeCoupledInertias, [("phi1", "phi2", "phi3"), ("w1","w2","w3"), ("tau0","tau1","tau2")])

end
