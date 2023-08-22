module TestMultiReturningFunction10

using Modia.StaticArrays
using Modia
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
        println("      TestMultiReturningFunction10: init Mbs with phi2=$phi2")
        ndof = length(phi2)
        obj = new(phi1, w1, 0.0, 0.0, phi2, zeros(ndof), zeros(ndof), zeros(ndof))
        # get_sizes(evaluatedParameters, ...)
        # add_states(obj,equationInfo, copy_from_x, copy_to_der_x)
        return obj
    end
end

struct Dummy
end;

function myBuildFunction(model::AbstractDict, modelModule, FloatType::Type, TimeType::Type,
                         instantiateModelOptions, ID, modelPathAST; buildOption = "Default")
    modelPathAsString = if isnothing(modelPathAST); "" else string(modelPathAST) end
    println("  TestMultiReturningFunction10: Test output from function myBuildFunction at modelPath = \"$modelPathAsString\":\n  Code could be constructed here and merged to the model with buildOption=$buildOption")
    return (model, Dummy())
end

MyModelWithBuild(; kwargs...) = Model(; _buildFunction = Par(functionName = :myBuildFunction), kwargs...)

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

Pendulum = MyModelWithBuild(_buildOption = "MyBuildOption",
    phi1 = Var(init=pi/2),
    w1   = Var(init=0.0),
    qdd  = Var(start=zeros(2)),
    mbs  = Mbs(phi2=[1.0,2.0]),
    mbs1 = Var(hideResult=true),
    mbs2 = Var(hideResult=true),
    mbs3 = Var(hideResult=true),
    mbs4 = Var(hideResult=true),
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
simulate!(pendulum, stopTime = 2.0, log=true)
showInfo(pendulum)

plot(pendulum, [("phi1", "w1"), "der(w1)", "qdd"])

simulate!(pendulum, stopTime = 2.0, log=true, merge=Map(mbs = Mbs(phi2=[10.0,20.0,30.0]), qdd = zeros(3)))
plot(pendulum, [("phi1", "w1"), "der(w1)", "qdd"], figure=2)

end
