module TestRectifier

println("\nRectifier: Demonstrating the ability to simulate models with state events")

using Modia
using Modia.Electric

# Desired:
#   using ModiaMath: plot
#
# In order that ModiaMath need not to be defined in the user environment, it is included via Modia:
using Modia.ModiaMath: plot


@model Rectifier begin
  R=Resistor(R=1)
  r=Resistor(R=0.01)
  C=Capacitor(C=1,start=0)
  D=IdealDiode()
  V=SineVoltage(V=5,freqHz=1.5, offset=0, startTime=0)
@equations begin
  connect(V.p, D.p)
  connect(D.n, R.p)
  connect(R.n, r.p)
  connect(r.n, V.n)
  connect(C.n, R.n)
  connect(C.p, R.p)
  end
end 

result = checkSimulation(Rectifier, 2, "C.v", 0.4773911315322196, logTranslation=true)

plot(result, ("C.v", "V.v"), heading="Rectifier", figure=12)

end
