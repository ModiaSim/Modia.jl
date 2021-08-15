module TestMultiReturningFunction5

using StaticArrays
using ModiaLang
@usingModiaPlot

mutable struct Mbs
    phi
    w
    derw
end

function setStates(phi, w)
    return Mbs(phi,w,0.0)
end

function setAccelerations(mbs::Mbs,derw)
    mbs.derw = derw
    return mbs
end

function getForces(mbs::Mbs)   
    # m*L^2*derw + m*g*l*sin(phi) ) = u
    #   m=L=g=1
    # derw + sin(phi) = 0
    u = mbs.derw + sin(mbs.phi)
    return SVector(u)
end


Pendulum = Model(
    d   = 1.0,   # damping constant
    phi = Var(init=pi/2),
    w   = Var(init=0.0),
    equations = :[
        w    = der(phi)
        mbs1 = setStates(phi,w)
        mbs2 = setAccelerations(mbs1,der(w))
        (u,) = getForces(mbs2)
        u = -d*w
    ]
)
    
pendulum = @instantiateModel(Pendulum , unitless=true, log=true, logDetails=false, logCode=true, logStateSelection=false)

simulate!(pendulum, stopTime = 2.0, log=true)

plot(pendulum, ("phi", "w"))

end
