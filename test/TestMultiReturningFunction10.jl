module TestMultiReturningFunction10

using ModiaLang.StaticArrays
using ModiaLang
@usingModiaPlot

mutable struct MbsData
    phi1::Float64
    w1::Float64
    derw1::Float64
    u1::Float64
    
    phi2::Vector{Float64}
    w2::Vector{Float64}
    derw2::Vector{Float64}
    residuals2::Vector{Float64}

    function MbsData(; path = "", phi1=0.0, w1=0.0, phi2, equationInfo = missing, evaluatedParameters=missing)
        println("... init Mbs with phi2=$phi2")
        ndof = length(phi2)
        obj = new(phi1, w1, 0.0, 0.0, phi2, zeros(ndof), zeros(ndof), zeros(ndof))
        # get_sizes(evaluatedParameters, ...)
        # add_states(obj,equationInfo, copy_from_x, copy_to_der_x)
        return obj
    end
end

Mbs(; kwargs...) = Par(; _constructor = :(MbsData), _path = true, kwargs...)

function copy_from_x(mbs::MbsData, x, ibeg, iend)
end

function copy_to_der_x(mbs::MbsData, der_x, ibeg, iend)
end

function setStates(mbs::MbsData, phi1, w1)::MbsData
    mbs.phi1 = phi1
    mbs.w1   = w1
    return mbs
end

function setAccelerations1(mbs::MbsData,derw1)::MbsData
    mbs.derw1 = derw1
    return mbs
end

function setAccelerations2(mbs::MbsData,derw2)::MbsData
    mbs.derw2 = derw2
    return mbs
end

function computeForcesAndResiduals(mbs::MbsData,time)
    # m*L^2*derw + m*g*l*sin(phi) ) = u
    #   m=L=g=1
    # derw + sin(phi) = 0
    mbs.u1         = mbs.derw1 + sin(mbs.phi1)
    mbs.residuals2 = (1+time)*mbs.derw2 + sin.(mbs.phi2)
    return mbs
end

function getForces(mbs::MbsData)
    return mbs.u1
end

function getResiduals(mbs::MbsData)
    return mbs.residuals2
end

Pendulum = Model(
    phi1 = Var(init=pi/2),
    w1   = Var(init=0.0),
    qdd  = Var(start=zeros(2)),
    mbs  = Mbs(phi2=[1.0,2.0]),
    equations = :[ 
        w1   = der(phi1)
        mbs1 = setStates(mbs,phi1,w1)
        mbs2 = setAccelerations1(mbs1,der(w1))
        mbs3 = setAccelerations2(mbs2,qdd)
        mbs4 = computeForcesAndResiduals(mbs3,time)
        tau1 = implicitDependency(getForces(mbs4), der(w1))
        0    = implicitDependency(getResiduals(mbs4), qdd)
        tau1 = -0.1*w1
    ]
)
        
pendulum = @instantiateModel(Pendulum , unitless=true, log=false, logDetails=false, logCode=true, logStateSelection=false)

#=
    function getDerivatives2(_der_x, _x, _m, _time)::Nothing
        _m.time = ModiaLang.getValue(_time)
        _m.nGetDerivatives += 1
        instantiatedModel = _m
        _p = _m.evaluatedParameters
        _leq_mode = nothing
        time = _time
        w1 = _x[1]
        phi1 = _x[2]
        mbs = Mbs()
        var"der(phi1)" = w1
        mbs1 = setStates(mbs, phi1, w1)
        tau1 = -0.1w1
        begin
            global var"der(w1)", qdd, mbs2, mbs3, mbs4
            _leq_mode = initLinearEquationsIteration!(_m, 1)
             ModiaBase.TimerOutputs.@timeit _m.timer "LinearEquationsIteration" while ModiaBase.LinearEquationsIteration!(_leq_mode, _m.isInitial, _m.solve_leq, _m.storeResult, _m.time, _m.timer)
                    var"der(w1)" = _leq_mode.x[1]
                    qdd = _leq_mode.x_vec[1]
                    mbs2 = setAccelerations1(mbs1, var"der(w1)")
                    mbs3 = setAccelerations2(mbs2, qdd)
                    mbs4 = computeForcesAndResiduals(mbs3, time)
                    ModiaBase.appendResidual!(_leq_mode.residuals, ustrip(getResiduals(mbs4)))
                    ModiaBase.appendResidual!(_leq_mode.residuals, ustrip(getForces(mbs4) - tau1))
                end
            _leq_mode = nothing
        end
        _der_x[1] = var"der(w1)"
        _der_x[2] = var"der(phi1)"
        if _m.storeResult     
            ModiaLang.addToResult!(_m, _der_x, time, mbs, mbs1, mbs2, mbs3, deepcopy(qdd), mbs4, tau1)
        end
        return nothing
    end

pendulum.getDerivatives! = getDerivatives2
=#


#simulate!(pendulum, stopTime = 0.1, interval = 0.05, log=true)
simulate!(pendulum, stopTime = 2.0, log=true)
# printResultInfo(pendulum)

plot(pendulum, [("phi1", "w1"), "der(w1)", "qdd"])


simulate!(pendulum, stopTime = 2.0, log=true, merge=Map(mbs = Mbs(phi2=[10.0,20.0,30.0]), qdd = zeros(3)))
plot(pendulum, [("phi1", "w1"), "der(w1)", "qdd"], figure=2)

end
