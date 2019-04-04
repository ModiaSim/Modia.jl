module SynchronousExamples

println("\nSynchronousExamples: Demonstrating the ability to simulate models with synchronous semantics")

using Modia
using Modia.Synchronous: sample, Clock, previous, hold

# Desired:
#   using ModiaMath: plot
#   using Test
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Modia.Test
end

@testset "Synchronous" begin

@model SynchronousOperators begin
  a = Var(size=())
  b = Var(size=())
  c = Var(start=0.0)
  obs=Float(start=0.0)  
@equations begin
  a = Clock(0.1)
  b = sample(time, a)
  c = previous(c, a) + 0.5
  der(obs) = 1000000*(c-obs)
  end
end

#=
result = simulate(SynchronousOperators, 5.0, storeEliminated=false, logSimulation=false)
plot(result, ("obs"), heading="SynchronousOperators",figure=15)
@test result["obs"][end] == 25
=#

@model MassWithSpringDamper begin
  m=1
  k=1 # spring constant
  d=0.1 # damping coefficient
  x=Float(start=0.0) # Position
  v=Float(start=0.0) # Velocity
  f=Float() # Force
@equations begin
  der(x) = v
  m*der(v) = f - k*x - d*v
  end
end 

@model SpeedControl begin
  @extends MassWithSpringDamper(k=0)
  @inherits v, f
  K = 5 # Gain of speed P controller
  vref = 100 # Speed ref.
  vd=Float(start=0.0)
  u=Float(start=0.0)
  fobs=Float(start=0.0)
@equations begin
  # speed sensor
  vd = sample(v, Clock(0.1))

  # P controller for speed
  u = K*(vref-vd)

  # force actuator
  f = hold(u)
  der(fobs) = 1000000*(f-fobs)
  end
end 

result = simulate(SpeedControl, 5.0, storeEliminated=false, logSimulation=true)
plot(result, ("v", "fobs"), heading="SpeedControl", figure=15)
@test result["v"][end] == 98.03921568627446

@model SpeedControlPI begin
  @extends MassWithSpringDamper(k=0)
  @inherits v, f
  K = 2 # Gain of speed P controller
  Ti = 2 # Integral time 
  vref = 100 # Speed reference
  dt=0.1 # sampling interval
  vd=Float()
  u=Float(start=0)
  e=Float()
  intE=Float(start=0)
    
  fobs=Float(start=0.0)
@equations begin
  # speed sensor
  vd = sample(v, Clock(dt))

  # PI controller for speed
  e = vref-vd
  intE = previous(intE, Clock(dt)) + e
  u = K*(e + intE/Ti)

  # force actuator
  f = hold(u)
  der(fobs) = 1000000*(f-fobs)
  end
end 

result = simulate(SpeedControlPI, 5.0, storeEliminated=false, logSimulation=false)
plot(result, ("v", "fobs"), heading="SpeedControlPI", figure=16)
@show result["v"][end]
@test isapprox(result["v"][end], 100.2849917097788; atol=1e-3)

@model ControlledMassBasic begin
  @extends MassWithSpringDamper(k=0) # k=100
  @inherits x, v, f
  KOuter = 10 # Gain of position PI controller
  KInner = 20 # Gain of speed P controller
  Ti = 10 # Integral time for pos. PI controller
  xref = 10 # Position reference
  xd=Float(size=())
  eOuter=Float(size=())
  intE=Float(start=0.0)
  uOuter=Float(size=())
  vd=Float(size=())
  vref=Float(size=())
  uInner=Float(start=0.0)
  fobs=Float(start=1990)
@equations begin
  # position sensor
  xd = sample(x, Clock(0.02))
  # outer PI controller for position
  eOuter = xref-xd
  intE = previous(intE, Clock(0.02)) + eOuter
  uOuter = KOuter*(eOuter + intE/Ti)
  # speed sensor
  vd = sample(v, Clock(0.02))
  # inner P controller for speed
  vref = uOuter
  uInner = KInner*(vref-vd)
  # force actuator
  f = hold(uInner)
  der(fobs) = 100000*(f-fobs)
  end
end

#result = simulate(ControlledMassBasic, 5.0, storeEliminated=false, logSimulation=true)
#plot(result, ("v", "fobs"), figure=4)


@model DiscretePIController begin
  K=0.1 # Gain 
  Ti=1E10 # Integral time 
  dt=1.0 # sampling interval
  ref=1 # set point
  u=Float(); ud=Float(size=())
  y=Float(); yd=Float(size=())
  e=Float(size=())
  intE=Float(start=0)

  fobs=Float(start=0)
@equations begin
  # sensor
  ud = sample(u, Clock(dt))
  # PI controller
  e = ref-ud
  intE = previous(intE, Clock(dt)) + e
  yd = K*(e + intE/Ti)
  # actuator
  y = hold(yd)
  
  der(fobs) = 100000*(y-fobs)
  end
end


@model SpeedControl2 begin
  m=MassWithSpringDamper(d=0.5, k=0)
  c=DiscretePIController(K=2, dt=0.1)
@equations begin
  connect(m.v, c.u)
  connect(c.y, m.f)
  end
end

#result = simulate(SpeedControl2, 5.0, storeEliminated=false, logSimulation=true)
#plot(result, ("v", "fobs"), figure=5)


@model SpeedControl3 begin
  m=MassWithSpringDamper(d=0.5, k=0)
  @extends DiscretePIController(K=2, dt=0.1)
  @inherits u, y, dt, e, intE
@equations begin
  u = m.v
  m.f = y
  end
end

#result = simulate(SpeedControl3, 5.0, storeEliminated=false, logSimulation=true)
#plot(result, ("v", "fobs"), figure=6)

end # testset
end
