

using Modia
using Base.Test


@testset "RunExamples" begin

Modia.ModiaLogging.resetTestStatus()

global figure = 1
include("CurrentController.jl")
include("Rectifier.jl")
@static if ! (VERSION < v"0.7.0-DEV.2005")
  include("CauerLowPassFilter.jl")
end
include("LinearSystems.jl")
include("SynchronousExamples.jl")
include("ElectricalVehicleAndCharger.jl")
include("CollidingBalls.jl")

Modia.ModiaLogging.printTestStatus()

end # testset
nothing

