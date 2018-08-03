

using Modia
using Base.Test


@testset "RunSimulations" begin

Modia.ModiaLogging.resetTestStatus()

include("models/TestVariableTypes.jl")
include("models/TestEquations.jl")

include("models/TestElectrical.jl")
include("models/TestFilter.jl")
include("models/TestArrayOfComponents.jl")

#include("models/TestCoupledInertias.jl")
#include("models/TestPendulum.jl")
#include("models/TestStateSelection.jl")
#include("models/TestRetranslationOmega.jl")

Modia.ModiaLogging.printTestStatus()

end # testset
nothing
