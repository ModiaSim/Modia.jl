

@static if VERSION < v"0.7.0-DEV.2005"
  using Base.Test
else
  using Test
end

@testset "RunTests" begin

include("runsymbolictests.jl")
#include("TestModelsWithError.jl")
include("runsimulationtests.jl")
include("../examples/runexamples.jl")
#include("runlargetests.jl")

end # testset
nothing
