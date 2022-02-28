module TestMultiReturningFunction7A

using Modia.StaticArrays
using Modia
@usingModiaPlot

mutable struct Mbs
    phi1
    w1
    phi2
    w2
    derw1
    derw2

    Mbs() = new(0,0,0,0,0,0)
end

function setStates(mbs, phi1, w1, phi2, w2)
    mbs.phi1 = phi1
    mbs.w1   = w1
    mbs.phi2 = phi2
    mbs.w2   = w2
    return mbs
end

function setAccelerations(mbs::Mbs,derw1,derw2)
    mbs.derw1 = derw1
    mbs.derw2 = derw2
    return mbs
end

function getForces(mbs::Mbs)
    # m*L^2*derw + m*g*l*sin(phi) ) = u
    #   m=L=g=1
    # derw + sin(phi) = 0
    u1 = mbs.derw1 + sin(mbs.phi1)
    u2 = mbs.derw2 + sin(mbs.phi2)
    return SVector(u1,u2)
end


Pendulum = Model(
    phi1 = Var(init=pi/2),
    w1   = Var(init=0.0),
    phi2 = Var(init=pi/4),
    w2   = Var(init=0.0),
    mbs  = Mbs(),
    equations = :[
        w1   = der(phi1)
        w2   = der(phi2)
        mbs1 = setStates(mbs,phi1,w1,phi2,w2)
        (0,0) = getForces(setAccelerations(mbs1,der(w1),der(w2)))
    ]
)

pendulum = @instantiateModel(Pendulum , unitless=true, log=false, logDetails=false, logCode=false, logStateSelection=false)

simulate!(pendulum, stopTime = 2.0, log=true,
          requiredFinalStates=[-1.399248030722412, -0.20564002545046062, -0.7148846477289476, -0.2742354977425878])

plot(pendulum, ("phi1", "w1", "phi2", "w2"))

end
