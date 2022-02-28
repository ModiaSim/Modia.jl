module MotorControl

println("\nMotorControl: Demonstrating the ability to simulate hierarchical mixed domain models")

using Modia
@usingModiaPlot
using  Modia.Measurements
export MotorControl2

include("$(Modia.modelsPath)/Blocks.jl")  
include("$(Modia.modelsPath)/Electric.jl")  
include("$(Modia.modelsPath)/Rotational.jl")  



MotorControl1 = Model(
    step = Step | Map(height=4.7*u"A", offset=0u"A"),
    feedback = Feedback,
    PI = PI | Map(T=0.005u"s", k=30, x = Var(init=0.0u"A")),
    firstOrder = FirstOrder | Map(k=1u"V/A", T=0.001u"s", x = Var(init=0.0u"V")),

    signalVoltage = SignalVoltage,
    resistor = Resistor | Map(R=13.8u"Ω"),
    inductor = Inductor | Map(L=0.061u"H"),
    emf = EMF | Map(k=1.016u"N*m/A"),
    ground = Ground,
    currentSensor = CurrentSensor,

    motorInertia = Inertia | Map(J=0.0025u"kg*m^2"),
    idealGear = IdealGear | Map(ratio=105),
    springDamper = SpringDamper | Map(c=5.0e5u"N*m/rad", d=500u"N*m*s/rad"),
    loadInertia = Inertia | Map(J=100u"kg*m^2"),
    tload = Torque,

    equations = :[
      tload.tau = 0.0u"N*m",
    ],

    connect = :[
      (step.y, feedback.u1)
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
      (motorInertia.flange_b, idealGear.flange_a)
      (idealGear.flange_b, springDamper.flange_a)
      (springDamper.flange_b, loadInertia.flange_a)
      (tload.flange, loadInertia.flange_b) ]
)

model = @instantiateModel(MotorControl1, log=false)
println("Simulate")
@time simulate!(model, stopTime=0.1, tolerance=1e-6, log=false, requiredFinalStates = 
    [3.487844601078223, 106.82874860310305, 4.616087152802267, 2.014120727821878, 41.98048886114646, 0.018332432934516536, 0.3725930536373392])
plot(model, [("currentSensor.i", "step.y"), "loadInertia.w"], figure=1)


# Hierarchical model

ControlledMotor = Model(
  refCurrent = input,
  feedback = Feedback,
  PI = PI | Map(T=0.005u"s", k=30, x = Var(init=0.0u"A")),
  firstOrder = FirstOrder | Map(k=1u"V/A", T=0.001u"s", x = Var(init=0.0u"V")),

  signalVoltage = SignalVoltage,
  resistor = Resistor | Map(R=13.8u"Ω"),
  inductor = Inductor | Map(L=0.061u"H"),
  emf = EMF | Map(k=1.016u"N*m/A"),
  ground = Ground,
  currentSensor = CurrentSensor,

  motorInertia = Inertia | Map(J=0.0025u"kg*m^2"),
  flange = Flange,

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

MotorControl2 = Model(
  step = Step | Map(height=4.7*u"A", offset=0u"A"),

  controlledMotor = ControlledMotor,

  idealGear = IdealGear | Map(ratio=105),
  springDamper = SpringDamper | Map(c=5.0e5u"N*m/rad", d=500u"N*m*s/rad"),
  loadInertia = Inertia | Map(J=100.0u"kg*m^2"),
  tload = Torque,

  equations = :[
    tload.tau = 0.0u"N*m",
  ],

  connect = :[
    (step.y, controlledMotor.refCurrent)

    (controlledMotor.flange, idealGear.flange_a)
    (idealGear.flange_b, springDamper.flange_a)
    (springDamper.flange_b, loadInertia.flange_a)
    (tload.flange, loadInertia.flange_b) ]
)

model = @instantiateModel(MotorControl2)
println("Simulate")
@time simulate!(model, stopTime=0.1, tolerance=1e-6, log=false, requiredFinalStates =  
    [3.487844601078223, 106.82874860310305, 4.616087152802267, 2.014120727821878, 41.98048886114646, 0.018332432934516536, 0.3725930536373392])
plot(model, [("controlledMotor.currentSensor.i", "step.y"), "loadInertia.w"], figure=1)


# Model with uncertainties

MotorControlWithUncertainties = MotorControl2 | Map( loadInertia = Map(J=(100.0 ± 10)*u"kg*m^2"), controlledMotor = Map(PI = Map(k=30 ± 3) ) )

model = @instantiateModel(MotorControlWithUncertainties, FloatType = Measurement{Float64})
println("Simulate")
@time simulate!(model, stopTime=0.1, tolerance=1e-6, log=false)
plot(model, [("controlledMotor.currentSensor.i", "step.y"), "loadInertia.w"], figure=1)

end

#=
module MotorControlModuleMonteCarlo

using Modia
using Unitful
using Main.MotorControlModule
include("../test/SimulateAndPlot.jl")

using MonteCarloMeasurements

# Model with Monte Carlo

MotorControlWithUncertainties = MotorControl2 | Map( loadInertia = Map(J=(100.0 ∓ 10)), controlledMotor = Map(PI = Map(k=30 ∓ 3) ) )

model = @instantiateModel(MotorControlWithUncertainties, FloatType = StaticParticles{Float64,100}, unitless=true)
println("Simulate")
@time simulate!(model, stopTime=0.1, tolerance=1e-6, log=false)
plot(model, [("controlledMotor.currentSensor.i", "step.y"), "loadInertia.w"], figure=1)

end
=#
