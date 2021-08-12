module TestMultiReturningFunction

using ModiaLang
@usingModiaPlot

include("$(ModiaLang.path)/models/Blocks.jl")
include("$(ModiaLang.path)/models/Electric.jl")
include("$(ModiaLang.path)/models/Rotational.jl")


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


    function getDerivatives(_der_x, _x, _m, _time)::Nothing
        _m.time = ModiaLang.getValue(_time)
        _m.nGetDerivatives += 1
        instantiatedModel = _m
        _p = _m.evaluatedParameters
        _leq_mode = nothing
        time = _time * upreferred(u"s")
        phi2 = _x[1]
        w2 = _x[2]
        tau = if time < 1 * u"s"
                _p[:tau_max]
            else
                if time < 2 * u"s"
                    0
                else
                    if time < 3 * u"s"
                        -(_p[:tau_max])
                    else
                        0
                    end
                end
            end
        phi1 = _p[:r] * phi2
        var"der(phi2)" = w2
        var"der(phi1)" = _p[:r] * var"der(phi2)"
        w1 = var"der(phi1)"
        begin
            local tau2, var"der(w2)", var"der(der(phi2))", var"der(der(phi1))", var"der(w1)", tau1
            _leq_mode = _m.linearEquations[1]
            _leq_mode.mode = -3
             ModiaBase.TimerOutputs.@timeit _m.timer "LinearEquationsIteration" while ModiaBase.LinearEquationsIteration(_leq_mode, _m.isInitial, _m.solve_leq, _m.storeResult, _m.time, _m.timer, useAppend = true)
                    tau2 = _leq_mode.x[1]
                    var"der(w2)" = tau2 / _p[:J2]
                    var"der(der(phi2))" = var"der(w2)"
                    var"der(der(phi1))" = _p[:r] * var"der(der(phi2))"
                    var"der(w1)" = var"der(der(phi1))"
                    tau1 = tau2 / _p[:r]
                    v = ustrip((tau - tau1) - _p[:J1] * var"der(w1)")
                    @show typeof(v)
                    append!(_leq_mode.residuals, v)
                end
            _leq_mode = nothing
        end
        _der_x[1] = ModiaLang.stripUnit(var"der(phi2)")
        _der_x[2] = ModiaLang.stripUnit(var"der(w2)")
        if _m.storeResult
            ModiaLang.addToResult!(_m, _der_x, time, tau, w1, var"der(phi1)", phi1, var"der(w1)", tau1, tau2, var"der(der(phi1))", var"der(der(phi2))")
        end
        return nothing
    end
end

threeCoupledInertias.getDerivatives!=getDerivatives

simulate!(threeCoupledInertias, stopTime = 2.0, log=true, 
          requiredFinalStates=[0.3029735305821086, 0.3933746028781303, 0.39337460287813036, 0.3029735305821086])

plot(threeCoupledInertias, [("phi1", "phi2", "phi3", "inertia1.phi"), 
                            ("w1","w2","w3", "inertia1.w"), 
                            ("tau0"),
                            ("tau1", "inertia2.flange_a.tau"),                            
                            ("tau2", "inertia3.flange_a.tau")])

end
