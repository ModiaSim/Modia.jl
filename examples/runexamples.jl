

#using Modia
@static if VERSION < v"0.7.0-DEV.2005"
  using Base.Test
else
  using Test
end

@testset "RunExamples" begin

Modia.ModiaLogging.resetTestStatus()

global figure = 1
include("CurrentController.jl")
include("Rectifier.jl")
@static if ! (VERSION < v"0.7.0-DEV.2005")
#  include("CauerLowPassFilter.jl")
end
include("LinearSystems.jl")
include("SynchronousExamples.jl")
@static if VERSION < v"0.7.0-DEV.2005"
  include("ElectricalVehicleAndCharger.jl") # Problem in 1.0
end
include("CollidingBalls.jl")
include("HeatTransfer2D.jl")

Modia.ModiaLogging.printTestStatus()

end # testset
nothing

