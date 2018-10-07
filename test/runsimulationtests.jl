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


@testset "RunSimulations" begin

    Modia.ModiaLogging.resetTestStatus()

    include("models/TestVariableTypes.jl")
    include("models/TestEquations.jl")

    include("models/TestElectrical.jl")
    include("models/TestFilter.jl")
    include("models/TestArrayOfComponents.jl")
    include("models/TestConditionalComponents.jl")
    include("models/TestConditionalEquations.jl")
    
    include("models/TestSpatialDiscretization.jl")

    include("models/MergingModifiers.jl")

    #include("models/TestCoupledInertias.jl")
    #include("models/TestPendulum.jl")
    #include("models/TestStateSelection.jl")
    #include("models/TestRetranslationOmega.jl")

    Modia.ModiaLogging.printTestStatus()

end 
nothing
