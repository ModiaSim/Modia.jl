module TestHeatTransfer2

using Modia
using Modia.Test
include("$(Modia.path)/models/HeatTransfer.jl")

HeatedRod2 = Model(
    fixedT     = FixedTemperature | Map(T = 493.15u"K"),
    fixedQflow = FixedHeatFlow,
    rod        = InsulatedRod2 | Map(L=1.0u"m", A=0.0004u"m^2", nT = 5), # 5 temperature nodes
    equations  = :[connect(fixedT.port, rod.port_a),
                   connect(rod.port_b , fixedQflow.port)]
)

# Model with unit
heatedRod2a = @instantiateModel(HeatedRod2, logCode=false, unitless=false)
simulate!(heatedRod2a, stopTime = 1e5, log=true, merge=Map(rod = Map(nT=8)), logParameters=true, logEvaluatedParameters=true)
showInfo(heatedRod2a)

@usingModiaPlot
plot(heatedRod2a, [("fixedT.port.T", "rod.T"), "rod.der(T)"], figure=1)

# Model without unit
heatedRod2b = @instantiateModel(HeatedRod2, logCode=false, unitless=true)
simulate!(heatedRod2b, stopTime = 1e5, log=true, merge=Map(rod = Map(nT=100)), logParameters=true, logEvaluatedParameters=true)
showInfo(heatedRod2b)
plot(heatedRod2b, [("fixedT.port.T", "rod.T"), "rod.der(T)"], figure=2)
end