module ServoSystemSimulation

println("\nServoSystem: Demonstrating the ability to simulate hierarchical mixed domain models")

using Modia
@usingModiaPlot

include("$(Modia.modelsPath)/Blocks.jl")
include("$(Modia.modelsPath)/Electric.jl")
include("$(Modia.modelsPath)/Rotational.jl")

setLogMerge(false)

Gear = Model(
    flange_a = Flange,
    flange_b = Flange,
    gear     = IdealGear_withSupport | Map(ratio = 105.0),
    fixed    = Fixed,
    spring   = Spring | Map(c=5.84e5u"N*m/rad"),
    damper1  = Damper | Map(d=500.0u"N*m*s/rad"),
    damper2  = Damper | Map(d=100.0u"N*m*s/rad"), 

    connect = :[
        (flange_a     , gear.flange_a)
        (fixed.flange , gear.support, damper2.flange_b)
        (gear.flange_b, damper2.flange_a, spring.flange_a, damper1.flange_a)
        (flange_b     , spring.flange_b, damper1.flange_b)]
)


ControlledMotor = Model(
  refCurrent = input,
  flange = Flange,
  feedback      = Feedback,
  PI            = PI | Map(k=30, T=1.0u"s"),
  firstOrder    = FirstOrder | Map(k=1.0, T=0.001u"s"),
  signalVoltage = UnitlessSignalVoltage,
  resistor      = Resistor | Map(R=13.8u"Ω"),
  inductor      = Inductor | Map(L=0.061u"H"),
  emf           = EMF | Map(k=1.016u"N*m/A"),
  ground        = Ground,
  currentSensor = UnitlessCurrentSensor,
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
    refSpeed = input, 
    motorSpeed = input,
    refCurrent = output,
    gain     = Gain | Map(k=105.0),
    PI       = PI | Map(T=1.0u"s", k=1.0),
    feedback = Feedback,

    connect = :[
        (refSpeed  , gain.u)
        (gain.y    , feedback.u1)
        (motorSpeed, feedback.u2)
        (feedback.y, PI.u)
        (PI.y      , refCurrent)]
)


Servo = Model(
    refSpeed        = input,
    flange_b        = Flange,
    speedController = SpeedController | Map(ks=1.0, Ts=1.0u"s", ratio=105.0),
    motor           = ControlledMotor | Map(km=30.0, Tm=0.005u"s"),
    gear            = Gear | Map(ratio=105.0),
    speedSensor1    = UnitlessSpeedSensor,
    speedSensor2    = UnitlessSpeedSensor,
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
    ks    = 0.4,
    Ts    = 0.08u"s",
    ramp  = Ramp  | Map(duration=1.18u"s", height=2.95),
    servo = Servo | Map(ks=:(up.ks), Ts=:(up.Ts)),
    load  = Inertia | Map(J=170.0u"kg*m^2"),
    equations =:[load.flange_b.tau = 0.0u"N*m"],
    connect = :[
        (ramp.y        , servo.refSpeed)
        (servo.flange_b, load.flange_a) ]
)

file = joinpath(pwd(), "TestServo.json")
writeModel(file, TestServo)
TestServo = readModel(file) | Map(ks = 0.8)

plotVariables = [("ramp.y", "load.w")         "servo.speedError.y";
                 "servo.gear.spring.phi_rel"  "servo.motor.currentSensor.i"]

testServo1 = @instantiateModel(TestServo, saveCodeOnFile="TestServo.jl")
println("Simulate")
@time simulate!(testServo1, Tsit5(), stopTime=2.0, tolerance=1e-6, requiredFinalStates = 
    [7.320842067204029, 9.346410309487013, 355.30389168655955, 2.792544498835712, 429.42665751348284, 311.7812493890421, 4.089776248793499, 2.969353608933471])
plot(testServo1, plotVariables, figure=1)

#=
println("\nServo with uncertainties")
using Modia.Measurements
TestServoWithUncertainties = TestServo | Map(
    load  = Map(J=(110.0 ± 20)u"kg*m^2"),
    servo = Map(motor = Map(resistor = Map(R=(13.8±1.0)u"Ω"))))
testServo2 = @instantiateModel(TestServoWithUncertainties , FloatType = Measurement{Float64})
@time simulate!(testServo2, Tsit5(), stopTime=2.0, tolerance=1e-6)
plot(testServo2, plotVariables, figure=2)

println("\nServo with Monte Carlo simulation")
using  Modia.MonteCarloMeasurements.Distributions
const nparticles = 100
uniform(vmin,vmax) = Modia.MonteCarloMeasurements.StaticParticles(nparticles,Distributions.Uniform(vmin,vmax))

TestServoWithMonteCarlo = TestServo | Map(
     load  = Map(J = uniform(50,170)),
     servo = Map(motor = Map(resistor = Map(R = uniform(12.8,14.8)))))
testServo3 = @instantiateModel(TestServoWithMonteCarlo, FloatType = Modia.MonteCarloMeasurements.StaticParticles{Float64,nparticles}, unitless=true)
@time simulate!(testServo3, Tsit5(), stopTime=2.0, tolerance=1e-6)
plot(testServo3, plotVariables, figure=3)
=#

end
