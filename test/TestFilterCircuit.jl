module TestFilterCircuit

using TinyModia

using DifferentialEquations
using ModiaPlot
using  Test

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

# Test access functions  
@testset "Test variable access functions" begin    
    @test ModiaPlot.getNames(filterCircuit) == ["C.C", "C.class", "C.i", "C.n.i", "C.n.v", "C.p.i", "C.p.v", "C.v",
"C.v", "R.R", "R.class", "R.i", "R.n.i", "R.n.v", "R.p.i", "R.p.v", "R.v", "V.V", "V.i", "V.n.i", "V.n.v", "V.p.i", "V.p.v", "V.v", "class", "der(C.v)", "ground.p.i", "ground.p.v", "time"]
    @test ModiaPlot.hasSignal(filterCircuit, "R.v")
    @test ModiaPlot.hasSignal(filterCircuit, "C.n.v")   
    @test ModiaPlot.hasSignal(filterCircuit, "R.R")
    @test ModiaPlot.hasSignal(filterCircuit, "ground.p.i")
    @test ModiaPlot.hasSignal(filterCircuit, "R.p.vv") == false
    @test isapprox(get_lastValue(filterCircuit, "R.v") , 2.5751560978893453u"V" )
    @test isapprox(get_lastValue(filterCircuit, "C.n.v"), 0.0u"V")
    @test isapprox(get_lastValue(filterCircuit, "R.R")  , 5.0u"Ω")
    @test isapprox(get_lastValue(filterCircuit, "ground.p.i"), 0.0)
end

# Test plotting of variables, zero variables, parameters 
plot(filterCircuit, [("R.v", "C.v"), ("R.R", "ground.p.i")])

end