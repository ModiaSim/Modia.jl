module ElectricalVehicleAndCharger

println("ElectricalVehicleAndCharger: Demonstrates the ability to change models from Julia.")

using Modia
using Modia.Electric

using Modia:@equation
using Modia:addEquation!
using Modia:deleteEquation!


# Desired:
#   using ModiaMath: plot
#
# In order that ModiaMath need not to be defined in the user environment, it is included via Modia:
using Modia.ModiaMath: plot


# -----------------------------------

function concatenateTimeSeries(result1, result2)
  result = result1
  for (name, values) in result
    result[name] = vcat(result1[name], result2[name])
  end
  return result
end

minute = 60
hour = 60*minute

@model Battery begin
  # Very simple model: Battery modelled as a capacitor.
  @extends Capacitor()
  @inherits C, v
  charge = Var()
@equations begin
  charge = C*v
  end
end

@model Charger begin
  # Very simple model: a constant voltage behind a resistor.
  v0=1000
  r=100
  V = ConstantVoltage(V=v0)
  R = Resistor(R=r)
  p = Pin()
  n = Pin()
  v = Var()
@equations begin
  connect(V.p, R.n)
  connect(V.n, n)
  connect(R.p, p)
  v = p.v - n.v
  end
end  
result = simulate(Charger, 1*hour)

@model ElectricVehicle begin
  # Very simple model: Battery modelled as a capacitor and the load as a resistor.
  B = Battery(C=1000)
  load = Resistor(R=1)
  p = Pin()
  n = Pin()
@equations begin
  connect(B.p, load.p)
  connect(B.n, load.n)
  connect(B.n, n)
  connect(B.p, p)
  end
end  
result = simulate(ElectricVehicle, 1*hour)
Bv = result["B.v"][end]  # Available voltage after driving used to intialize next model

Bv = 10.0

@model ElectricalVehicleWithCharger begin
  ev = ElectricVehicle(B = Battery(C=1000, v = Float(start=Bv)))
  c = if rand(1:2)==1; (println("Super Charger"); Charger()) else (println("Standard Charger"); Charger(v0=500, r=200)) end 
 # y = Var(size=())
end

global now = 0
then = 1*hour

global result = simulate(ElectricalVehicleWithCharger, then)
now = then
Bv = result["ev.B.v"][end]

for i in 1:5
  global now
  global result
  println("\nConnecting Electric Vehicle to Charger")   
  #addEquation!(ElectricalVehicleWithCharger, :(connect(ev.p, c.p)))
  #addEquation!(ElectricalVehicleWithCharger, :(connect(ev.n, c.n)))
  addEquation!(ElectricalVehicleWithCharger, @equation(connect(this.ev.p, this.c.p)))
  addEquation!(ElectricalVehicleWithCharger, @equation(connect(this.ev.n, this.c.n)))
  #addEquation!(ElectricalVehicleWithCharger, @equation(this.y = this.ev.C.v))

  then = now+15*minute
  resultCharging = simulate(ElectricalVehicleWithCharger, then, startTime=now, logStatistics=false)
  now = then

  result = concatenateTimeSeries(result, resultCharging)
  plot(result, ("ev.B.charge"), figure=17)
  Bv = result["ev.B.v"][end]

  println("\nDisconnect")
  deleteEquation!(ElectricalVehicleWithCharger, @equation(connect(this.ev.p, this.c.p)))
  deleteEquation!(ElectricalVehicleWithCharger, @equation(connect(this.ev.n, this.c.n)))

  then = now+1*hour
  resultDriving = simulate(ElectricalVehicleWithCharger, then, startTime=now, logStatistics=false)
  now = then

  result = concatenateTimeSeries(result, resultDriving)
  plot(result, ("ev.B.charge"), figure=17)
  Bv = result["ev.B.v"][end]
end


end
