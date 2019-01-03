
# Desired:
#   using Test
#
# In order that Test need not to be defined in the user environment, it is included via Modia:
import Modia

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Modia.Test
end

@testset "RunTests" begin

    include("runsymbolictests.jl")
    # include("TestModelsWithError.jl")
    include("runsimulationtests.jl")
    include("../examples/runexamples.jl")
    # include("runlargetests.jl")

end # testset
nothing
