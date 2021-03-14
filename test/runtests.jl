module Runtests

using Test

@testset "Test TinyModia with simulation" begin

    include("TestVariables.jl") 
    include("TestFirstOrder.jl")    
    include("TestStateSelection.jl") 
    
    include("TestUnits.jl")  
    include("TestUncertainties.jl")
    include("TestUnitsAndUncertainties.jl")
    
    include("TestTwoInertiasAndIdealGear.jl")   
    include("TestTwoInertiasAndIdealGearWithUnits.jl")   
    include("TestTwoInertiasAndIdealGearWithUnitsAndUncertainties.jl")   
    include("TestTwoInertiasAndIdealGearWithMonteCarlo.jl")   
    #include("TestTwoInertiasAndIdealGearWithUnitsAndMonteCarlo.jl")  # MonteCarlo and Unitful not yet supported
    
    include("TestSingularLRRL.jl")  

    include("TestStateSpace.jl")  
    include("TestHeatTransfer.jl")  
    
    include("../examples/runexamples.jl")    

end

end