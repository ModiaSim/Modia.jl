module ConditionalComponent

println("\nRectifier: Demonstrating conditional components")

using Modia
using Modia.Electric
using ModiaMath: plot

@model Conditional begin
  R=Resistor(R=1)
  r=Resistor(R=0.01)
  C=Capacitor(C=1,start=0)
  D=IdealDiode()
  V=SineVoltage(V=5,freqHz=1.5, offset=0, startTime=0)
  inc = true
  R2=if inc; Resistor(R=1) end
@equations begin
  connect(V.p, D.p)
  connect(D.n, R.p)
  connect(R.n, r.p)
  connect(r.n, V.n)
  connect(C.n, R.n)
  connect(C.p, R.p)
  connect(C.p, R2.p)
  connect(C.n, R2.n)
  end
end 

result = simulate(Conditional, 2, logTranslation=true)

plot(result, ("C.v", "V.v"), heading="Conditional", figure=12)


end
