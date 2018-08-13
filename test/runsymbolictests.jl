#= Run as:
cd(".../Modia/test")
include("runsymbolictests.jl")
=#

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

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

