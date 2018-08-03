#= Run as:
cd(".../Modia/test")
include("runsymbolictests.jl")
=#

using Base.Test

@testset "SymbolicTests" begin

include("symbolic/BLTandPantelides/setup.jl")
include("symbolic/BLTandPantelides/testBLTandPantelides.jl")

include("symbolic/DAEquations/setup.jl")
include("symbolic/DAEquations/testSymbolicTransform.jl")

include("symbolic/TestExactlyRemoveSingularities.jl")

include("symbolic/TestStateSelectionAlgorithm.jl")

include("symbolic/TestTearingAlgorithm.jl")


end
nothing

