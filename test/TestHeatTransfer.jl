module TestHeatTransfer

using TinyModia
using Unitful
@usingModiaPlot
using DifferentialEquations

include("$(TinyModia.path)/models/HeatTransfer.jl")

SimpleRod = Model(
    fixedT = FixedTemperature | Map(T = 393.15u"K"),
         C = HeatCapacitor    | Map(C = 200u"J/K"),
         G = ThermalConductor | Map(G = 100u"W/K"),
    connect = :[
        (fixedT.port, G.port_a),
        (G.port_b   , C.port)
    ]
)

simpleRod = @instantiateModel(SimpleRod)
simulate!(simpleRod, Tsit5(), stopTime = 15, requiredFinalStates = [393.09468783194814])
plot(simpleRod, ("fixedT.port.T", "C.T"), figure=1)


HeatedRod = Model(
    fixedT     = FixedTemperature | Map(T = 493.15u"K"),
    fixedQflow = FixedHeatFlow, 
    rod = InsulatedRod | Map(T = Var(init = fill(293.15u"K", 5))),
    connect = :[
        (fixedT.port, rod.port_a),
        (rod.port_b , fixedQflow.port)
    ]
)

heatedRod = @instantiateModel(HeatedRod)
simulate!(heatedRod, Tsit5(), stopTime = 1e5, 
          requiredFinalStates = [492.9629728925529, 492.607421697723, 492.30480548105595, 492.0850627579106, 491.9694736486912])
plot(heatedRod, ("fixedT.port.T", "rod.T"), figure=2)
    
end
