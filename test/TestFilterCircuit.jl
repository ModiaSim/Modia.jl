module TestFilterCircuit

using Modia
@usingModiaPlot
using Modia.Test

setLogMerge(false)

include("../models/Electric.jl")

FilterCircuit = Model(
    R = Resistor | Map(R=0.5u"Ω"),
    C = Capacitor | Map(C=2.0u"F", v=Var(init=0.1u"V")),
    V = ConstantVoltage | Map(V=10.0u"V"),
    ground = Ground,
    connect = :[
      (V.p, R.p)
      (R.n, C.p)
      (C.n, V.n, ground.p)
    ]
)

filterCircuit = @instantiateModel(FilterCircuit)

simulate!(filterCircuit, Tsit5(), stopTime = 10, merge = Map(R = Map(R = 5u"Ω"), C = Map(v = 3.0u"V")), 
          logParameters = true, logStates = true, requiredFinalStates = [7.424843902110655]) 
showInfo(filterCircuit)

# Test access functions  
@testset "Test variable access functions (TestFilterCircuit.jl)" begin  
    currentNames  = getSignalNames(filterCircuit, getPar=false, getMap=false)
    requiredNames = String["C.i", "C.n.i", "C.n.v", "C.p.i", "C.p.v", "C.v", "R.i", "R.n.i", "R.n.v", "R.p.i", "R.p.v", "R.v", "V.i", "V.n.i", "V.n.v", "V.p.i", "V.p.v", "V.v", "C.der(v)", "ground.p.i", "ground.p.v", "time"]
    @test sort!(currentNames) == sort!(requiredNames)

    @test hasSignal(filterCircuit, "R.v")
    @test hasSignal(filterCircuit, "C.n.v")   
    @test hasParameter(filterCircuit, "R.R")
    @test hasSignal(filterCircuit, "ground.p.i")
    @test hasSignal(filterCircuit, "R.p.vv") == false
    @test isapprox(getLastValue(filterCircuit, "R.v") , 2.5751560978893453u"V" )
    @test isapprox(getLastValue(filterCircuit, "C.n.v"), 0.0u"V")
    @test isapprox(getEvaluatedParameter(filterCircuit, "R.R")  , 5.0u"Ω")
    @test isapprox(getLastValue(filterCircuit, "ground.p.i"), 0.0)
end

# Test plotting of variables, zero variables, parameters 
plot(filterCircuit, [("R.v", "C.v"), "ground.p.i"], figure=1)


# Simulate with lower precision
filterCircuitLow = @instantiateModel(FilterCircuit, unitless=true, FloatType = Float32)  # unitless=true required, since Float32 not yet supported, if units are present
simulate!(filterCircuitLow, RK4(), adaptive=false, stopTime=10.0, interval=0.01, 
          merge = Map(R = Map(R = 5.0), C = Map(v = 3.0)),                              #    merge = Map(R = Map(R = 5u"Ω"), C = Map(v = 3.0u"V")),
          requiredFinalStates = Float32[7.4248414])
plot(filterCircuitLow, [("R.v", "C.v"), "ground.p.i"], figure=2)


# Simulate with DAE integrator
println("\n... Simulate with DAE integrator")
filterCircuit2 = @instantiateModel(FilterCircuit) # Without a new @instantiateModel, this fails with a type error (will be fixed at a later time instant)
simulate!(filterCircuit2, IDA(), stopTime = 10, merge = Map(R = Map(R = 5u"Ω"), C = Map(v = 3.0u"V")), requiredFinalStates = [7.424843902110655]) 
    
end