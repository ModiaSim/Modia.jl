module TestSynchronous

#println("\nSynchronousExamples: Demonstrating the ability to simulate models with synchronous semantics")

using TinyModia, ModiaPlot
using DifferentialEquations


MassWithSpringDamper = Model(
    m = 1.0,
    k = 1.0,           # spring constant
    d = 0.1,           # damping coefficient
    x = Var(init=0.0), # Position
    v = Var(init=0.0), # Velocity
    equations1 = :[
        der(x) = v
        m*der(v) = f - k*x - d*v
    ]
) 

SpeedControl = MassWithSpringDamper | Model(
    k    = 0.0,
    K    = 5.0,   # Gain of speed P controller
    vref = 100.0, # Speed ref.
    vd = Var(start=0.0),
    u  = Var(start=0.0),  
    equations2 = :[
        c  = Clock(0.1, instantiatedModel, 1)
        vd = sample(v, c, instantiatedModel, 1)    # speed sensor
        u  = K*(vref-vd)                           # P controller for speed
        f  = hold(u)                               # force actuator
    ]
)

speedControl = @instantiateModel(SpeedControl)
simulate!(speedControl, Tsit5(), stopTime=1.5, log=true, logEvents=true)
plot(speedControl, [("v", "f"), "vd"], heading="SpeedControl", figure=1)

end