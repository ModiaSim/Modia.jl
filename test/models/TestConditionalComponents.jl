module ConditionalComponent

println("\nRectifier: Demonstrating conditional components")

using Modia
using Modia.Electric

# Desired:
#   using ModiaMath: plot
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot

@model ConditionalLoad begin
    R=Resistor(R=1)
    r=Resistor(R=0.01)
    C=Capacitor(C=1,start=0)
    D=IdealDiode()
    V=SineVoltage(V=5,freqHz=1.5, offset=0, startTime=0)
    extraLoad = true
    R2=if extraLoad; Resistor(R=1) end
@equations begin
    connect(V.p, D.p)
    connect(D.n, R.p)
    connect(R.n, r.p)
    connect(r.n, V.n)
    connect(C.n, R.n)
    connect(C.p, R.p)
    connect(C.n, R2.n)
    connect(C.p, R2.p)
    end
end 

result = simulate(ConditionalLoad, 2, logTranslation=true)
plot(result, ("C.v", "V.v"), heading="Conditional", figure=12)


@model NoExtraLoad begin
    @extends ConditionalLoad(extraLoad = false)
end

result = simulate(NoExtraLoad, 2) # , logTranslation=true) # Problem with logging nothing in Julia 1.0
plot(result, ("C.v", "V.v"), heading="Conditional", figure=12)

end
