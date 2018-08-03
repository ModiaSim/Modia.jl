

using Base.Test

@testset "RunTests" begin

include("runsymbolictests.jl")
#include("TestModelsWithError.jl")
include("runsimulationtests.jl")
include("../examples/runexamples.jl")
#include("runlargetests.jl")

end # testset
nothing
