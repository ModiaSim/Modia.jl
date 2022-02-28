module TestMultiReturningFunction5A

using Modia
@usingModiaPlot

mutable struct Mbs
    phi
    w
    derw
end

function setStates(phi, w)
    return Mbs(phi,w,0.0)
end

function setAccelerations(mbs::Mbs,derw)::Mbs
    mbs.derw = derw
    return mbs
end

function getForces(mbs::Mbs)   
    # m*L^2*derw + m*g*l*sin(phi) ) = u
    #   m=L=g=1
    # derw + sin(phi) = 0
    u = mbs.derw + sin(mbs.phi)
    return u
end


Pendulum = Model(
    d   = 1.0,   # damping constant
    phi = Var(init=pi/2),
    w   = Var(init=0.0),
    equations = :[
        w    = der(phi)
        mbs1 = setStates(phi,w)
        mbs2 = setAccelerations(mbs1,der(w))
        u    = implicitDependency(getForces(mbs2),der(w)) 
        u    = -d*w
    ]
)
    
pendulum = @instantiateModel(Pendulum , unitless=true, log=false, logDetails=false, logCode=false, logStateSelection=true)

simulate!(pendulum, stopTime = 10.0, log=true, requiredFinalStates= [-0.012169911296941314, 0.0024087599260249094])

plot(pendulum, ("phi", "w"))

end
