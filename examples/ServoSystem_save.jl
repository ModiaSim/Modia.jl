module ServoSystemSimulation

println("\nServoSystem: Demonstrating the ability to simulate hierarchical mixed domain models")

include("../models/Blocks.jl")  
include("../models/Electric.jl")  
include("../models/Rotational.jl")  

using TinyModia
using DifferentialEquations
using ModiaPlot
using Unitful
#using Measurements
using MonteCarloMeasurements


setLogMerge(false)

Gear = Model(
    flange_a = Flange,
    flange_b = Flange,
    gear     = IdealGear_withSupport | Map(ratio = 105),
    fixed    = Fixed,
    spring   = Spring | Map(c=5.84e5u"N*m/rad"),
    damper1  = Damper | Map(d=500u"N*m*s/rad"),
    damper2  = Damper | Map(d=100u"N*m*s/rad"),

    connect = :[
        (flange_a     , gear.flange_a)
        (fixed.flange , gear.support, damper2.flange_b)
        (gear.flange_b, damper2.flange_a, spring.flange_a, damper1.flange_a)
        (flange_b     , spring.flange_b, damper1.flange_b)]
)


ControlledMotor = Model(
  inputs = :[refCurrent],
  flange = Flange,
  feedback      = Feedback,
  PI            = PI | Map(k=30, T=1.0u"s", init = Map(x=0.0u"A")),
  firstOrder    = FirstOrder | Map(k=1u"V/A", T=0.001u"s", init = Map(x=0.0u"V")),
  signalVoltage = SignalVoltage,
  resistor      = Resistor | Map(R=13.8u"Ω"),
  inductor      = Inductor | Map(L=0.061u"H"),
  emf           = EMF | Map(k=1.016u"N*m/A"),
  ground        = Ground,
  currentSensor = CurrentSensor,
  motorInertia  = Inertia | Map(J=0.0025u"kg*m^2"),

  connect = :[
    (refCurrent, feedback.u1)
    (feedback.y, PI.u)
    (PI.y, firstOrder.u)
    (firstOrder.y, signalVoltage.v)

    (signalVoltage.p, resistor.p)
    (resistor.n, inductor.p)
    (inductor.n, emf.p)
    (emf.n, ground.p, currentSensor.p)
    (currentSensor.n, signalVoltage.n)
    (currentSensor.i, feedback.u2)

    (emf.flange, motorInertia.flange_a)
    (motorInertia.flange_b, flange) ]
)


SpeedController = Model(
    inputs   = :[refSpeed, motorSpeed],
    outputs  = :[refCurrent],
    gain     = Gain | Map(k=105),
    PI       = PI | Map(T=1.0u"s", k=1.0u"A*s/rad", init = Map(x=0u"rad/s")),
    feedback = Feedback,

    connect = :[
        (refSpeed  , gain.u)
        (gain.y    , feedback.u1)
        (motorSpeed, feedback.u2)
        (feedback.y, PI.u)
        (PI.y      , refCurrent)]
)


Servo = Model(
    inputs          = :[refSpeed],
    flange_b        = Flange,
    speedController = SpeedController | Map(ks=1.0u"A*s/rad", Ts=1.0u"s", ratio=105),
    motor           = ControlledMotor | Map(km=30.0, Tm=0.005u"s"),
    gear            = Gear | Map(ratio=105),
    speedSensor1    = SpeedSensor,
    speedSensor2    = SpeedSensor,
    speedError      = Feedback,
    
    connect = :[
        (refSpeed                  , speedController.refSpeed, speedError.u1)
        (speedController.refCurrent, motor.refCurrent)
        (motor.flange              , gear.flange_a, speedSensor1.flange)
        (speedSensor1.w            , speedController.motorSpeed)
        (gear.flange_b             , flange_b, speedSensor2.flange)
        (speedSensor2.w            , speedError.u2)]
)


TestServo = Model(
    ks    = 0.8u"A*s/rad",
    Ts    = 0.08u"s",
    ramp  = Ramp | Map(duration=1.18u"s", height=2.95u"rad/s", offset=0.0u"rad/s"),
    servo = Servo | Map(ks=:(up.ks), Ts=:(up.Ts)),
    load  = Inertia | Map(J=170u"kg*m^2"),
    equations =:[load.flange_b.tau = 0u"N*m"],
    connect = :[
        (servo.refSpeed, ramp.y)
        (load.flange_a, servo.flange_b) ]
)
    
plotVariables = [("ramp.y", "load.w"), "servo.speedError.y", "servo.motor.currentSensor.i"]

model = @instantiateModel(TestServo, log=false, logDetails=false)
println("Simulate")
@time simulate!(model, Tsit5(), stopTime=2.0, tolerance=1e-6, log=false)
plot(model, plotVariables, figure=1)


println("\nModel with uncertainties")
#TestServoWithUncertainties = TestServo | Map( load = Map(J=(110.0 ± 20)*u"kg*m^2") )
TestServoWithUncertainties = TestServo | Map( load = Map(J=(110.0 ∓ 20) ) ) # u"kg*m^2") ) #Map( ramp = Map(height=3 ∓ 1))

#model = @instantiateModel(TestServoWithUncertainties , FloatType = Measurement{Float64})
model = @instantiateModel(TestServoWithUncertainties , FloatType = StaticParticles{Float64,100}, unitless=true, log=false, 
    logCode=false, logExecution=false)
println("Simulate")
@time simulate!(model, Tsit5(), stopTime=2.0, tolerance=1e-6, log=false)
plot(model, plotVariables, figure=2)


end
