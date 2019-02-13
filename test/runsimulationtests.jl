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
    @static if VERSION < v"0.7.0-DEV.2005"
        # Don't run TestFilter since aliasElimination=true does not work under 0.6
    else
        include("models/TestFilter.jl")
    end
    include("models/TestArrayOfComponents.jl")
    include("models/TestConditionalComponents.jl")
    include("models/TestConditionalEquations.jl")
    
    include("models/TestSpatialDiscretization.jl")

    include("models/MergingModifiers.jl")

    #include("models/TestFluid.jl")

    #include("models/TestCoupledInertias.jl")
    #include("models/TestPendulum.jl")
    include("models/TestTearing.jl")
    include("models/TestAutomaticStateSelection.jl")
    #include("models/TestRetranslationOmega.jl")

    Modia.ModiaLogging.printTestStatus()

end 
nothing
