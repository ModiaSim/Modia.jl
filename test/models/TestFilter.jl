module TestFilter

println("TestFilter: Tests various features of the symbolic handling.")

using Modia
using Modia.Electric

# Desired:
#   using ModiaMath: plot
#
# In order that these packages need not to be defined in the user environment, they are included via Modia:
using Modia.ModiaMath: plot


@model LPfilter begin
    R = Resistor(R=100.0)
    C = Capacitor(C=2.5E-3, v=Float(start=0.0))
    V = ConstantVoltage(V=10.0)
    ground = Ground()
    @equations begin
        connect(V.n, ground.p)
        connect(V.p, R.p)
        connect(R.n, C.p)
        connect(C.n, V.n)
    end
end 
result = simulate(LPfilter, 2)
plot(result, ("R.v", "C.v", "V.v"), figure=3)

checkSimulation(LPfilter, 2, "C.v", 9.996843043981416, logTranslation=true, logSimulation=true, storeEliminated=false)
simulate(LPfilter, 2, aliasElimination=true)
checkSimulation(LPfilter, 2, "C.v", 9.996843043981416, aliasElimination=true, logName="LPfilter aliasElimination")
checkSimulation(LPfilter, 2, "C.v", 9.996843043981416, aliasElimination=true, removeSingularities=true, logName="LPfilter aliasElimination removeSingularities")
simulate(LPfilter, 2.0, aliasElimination=true, removeSingularities=true, logName="LPfilter aliasElimination removeSingularities")
#checkSimulation(LPfilter, 2, "C.v", 9.996843043981416, aliasElimination=true, removeSingularities=true, expandArrayIncidence=true, logName="LPfilter aliasElimination removeSingularities expandArrayIncidence")
#checkSimulation(LPfilter, 2, "C.v", 9.996843043981416, aliasElimination=true, removeSingularities=true, expandArrayIncidence=true, useIncidenceMatrix=true, logName="LPfilter aliasElimination removeSingularities expandArrayIncidence useIncidenceMatrix")
#checkSimulation(LPfilter, 2, "C.v", 9.996843043981416, newStateSelection=true)


@model LPfilterWithoutGround begin
    R = Resistor(R=100.0)
    C = Capacitor(C=2.5E-3)
    V = ConstantVoltage(V=10.0)
    @equations begin
        connect(V.p, R.p)
        connect(R.n, C.p)
        connect(C.n, V.n)
    end
end 

#res = simulate(LPfilterWithoutGround , 2)
checkSimulation(LPfilterWithoutGround, 2, "C.v", 9.996843043981416, removeSingularities=true)

@model LPfilterWithoutGround1 begin
    R = Resistor(R=100.0Ohm) # , n=Pin(v=Float(start=0.0)))
    C = Capacitor(C=2.5Milli * Farad, v=Voltage(start=0.0Volt))
    V = ConstantVoltage(V=10.0Volt)
    @equations begin
        connect(V.p, R.p)
        connect(R.n, C.p)
        connect(C.n, V.n)
    end
end 
#checkSimulation(LPfilterWithoutGround1, 2, "C.v", 9.996843043981416)

@model LPfilterWithoutGround2 begin
    R = Resistor(R=100) # , n=Pin(v=Float(start=0.0)))
    C = Capacitor(C=2.5E-3, v=Float(start=0.0))
    V = ConstantVoltage(V=10)
    @equations begin
        connect(V.p, R.p)
        connect(R.n, C.p)
        connect(C.n, V.n)
    end
end 


@model LPfilterAndSineSource begin
    R = Resistor(R=1)
    C = Capacitor(C=1)
    V = SineVoltage(V=1)
    @equations begin
        connect(V.p, R.p)
        connect(R.n, C.p)
        connect(C.n, V.n)
    end
end 
result = simulate(LPfilterAndSineSource, 10)
plot(result, ("C.v", "V.v"), heading="LPfilterAndSineSource", figure=4)


@model LPfilterModel begin
    p = Pin()
    n = Pin()
    R = Resistor(R=1)
    C = Capacitor(C=1)
    ground = Ground()
    @equations begin
        connect(p, R.p)
        connect(n, R.n)
        connect(R.n, C.p)
        connect(C.n, ground.p)
    end
end 

# Redeclaration of component models
@model HPfilter begin
    @extends LPfilter(R=Capacitor(C=2), C=Resistor(R=3))
end
checkSimulation(HPfilter, 2, "C.v", 7.16540372163548, removeSingularities=true)


# Redeclaration of models
@model GenericFilter begin
    m = Resistor
    @extends LPfilter(C=m(C=2E-3))
end

@model NewFilter begin
    @extends GenericFilter(m=Capacitor)
end
checkSimulation(NewFilter, 2, "C.v", 9.999596486928553, removeSingularities=true)
 
 
#high = true
# Conditional model constructor
@model CondFilter begin
    high = true  # Does not work without this.high.
    @extends LPfilter(R=if high Capacitor(C=2) else Resistor(R=3) end, 
    C=if high Resistor(R=1) else Capacitor(C=1) end)
end
checkSimulation(CondFilter, 2, "C.v", 3.678778096777374, removeSingularities=true)


# Conditional model
@model CondFilter2 begin
    high = true 
    @extends LPfilter(R=high ? Capacitor(C=2) : Resistor(R=1), 
    C=high ? Resistor(R=1) : Capacitor(C=2))
end
checkSimulation(CondFilter2, 2, "C.v", 3.678778096777374, removeSingularities=true)


# Model arrays
Models = [Resistor, Capacitor]
high = true
@model FilterModels begin
    selectModel = 2
    @extends LPfilter(R=high ? Models[selectModel](C=2) : Resistor(R=2), C=high ? Resistor(R=1) : Capacitor(C=2))
    @equations begin
  end
end
checkSimulation(FilterModels, 2, "C.v", 3.678778096777374, removeSingularities=true)


# Model component arrays 
ModelComponents = [Resistor(R=5), Capacitor(C=1)]

@model FilterComponents begin
    selectModel = 2
    high = true
    @extends LPfilter(R=high ? ModelComponents[selectModel] : Resistor(), 
    C=high ? Resistor(R=2) : Capacitor())
end
checkSimulation(FilterComponents, 2, "C.v", 3.678778096777374)


@model TenCoupledFilters begin
    V = ConstantVoltage(V=1)
    ground = Ground()
    F1 = LPfilterModel()
    F2 = LPfilterModel()
    F3 = LPfilterModel()
    F4 = LPfilterModel()
    F5 = LPfilterModel()
    F6 = LPfilterModel()
    F7 = LPfilterModel()
    F8 = LPfilterModel()
    F9 = LPfilterModel()
    F10 = LPfilterModel()
    @equations begin
        connect(V.n, ground.p)
        connect(V.p, F1.p)
        connect(F1.n, F2.p)
        connect(F2.n, F3.p)
        connect(F3.n, F4.p)
        connect(F4.n, F5.p)
        connect(F5.n, F6.p)
        connect(F6.n, F7.p)
        connect(F7.n, F8.p)
        connect(F8.n, F9.p)
        connect(F9.n, F10.p)
    end
end
checkSimulation(TenCoupledFilters, 2, "F10.C.v", 1.232726022885833e-5)
checkSimulation(TenCoupledFilters, 2, "F10.C.v", 1.232726022885833e-5, aliasElimination=true)

end
