module TestCurrentController

println("\nCurrentController: Demonstrating the ability to simulate mixed domain models")

using Modia
using Modia.Electric
using Modia.Rotational
using Modia.Blocks

# Desired:
#   using ModiaMath: plot
#
# In order that ModiaMath need not to be defined in the user environment, it is included via Modia:
using Modia.ModiaMath: plot


@model CurrentController begin
  k=30   # Gain of PI current controller
  T=0.005   # Time constant of PI current controller
  Jmotor=Inertia(J=0.0025) 
  emf=EMF(k=1.016) #, phi=Float(state=false))
  inductor=Inductor(L=0.061) 
  resistor=Resistor(R=13.8) 
  ground=Ground()
  signalVoltage=SignalVoltage()
  currentSensor=CurrentSensor()
  firstOrder=FirstOrder(k=1,T=0.001)
  feedback=Feedback()
  PI1=PI(T=T, k=k)
  gear=IdealGear(ratio=105)
  load=Inertia(J=100)
  step=Step(height=1,startTime=0)
  spring=SpringDamper(c=5e5, d=500)
@equations begin
  connect(emf.flange, Jmotor.flange_a)
  connect(resistor.n, inductor.p) 
  connect(ground.p, emf.n)
  connect(inductor.n, emf.p)
  connect(signalVoltage.p, resistor.p)
  connect(currentSensor.p, ground.p)
  connect(currentSensor.n, signalVoltage.n)
  connect(currentSensor.i, feedback.u2) 
  connect(PI1.y, firstOrder.u)
  connect(feedback.y, PI1.u) 
  connect(signalVoltage.v, firstOrder.y)
  connect(feedback.u1, step.y)
  connect(Jmotor.flange_b, gear.flange_a)
  connect(load.flange_a, spring.flange_b)
  connect(gear.flange_b, spring.flange_a)
  end
end

# result = simulate(CurrentController, 0.1, tearing=true, logTranslation=true) # Fails
result = checkSimulation(CurrentController, 0.1, "load.w", 0.07929151315487117, removeSingularities=false, tearing=true)

result = checkSimulation(CurrentController, 0.1, "load.w", 0.07927285979038304)
plot(result, [("currentSensor.i", "step.y"), "load.w"], heading="CurrentController", figure=11)

end
