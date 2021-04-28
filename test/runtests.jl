module Runtests

using Test

@testset "Test TinyModia with simulation" begin

    include("TestVariables.jl") 
    include("TestFirstOrder.jl") 
    include("TestFirstOrder2.jl")        
    include("TestStateSelection.jl") 
    include("TestFilterCircuit.jl")     
    
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
    include("TestParameter.jl") 
    include("TestHeatTransfer.jl")  
 
    include("TestSimpleStateEvents.jl") 
    include("TestSynchronous.jl") 
     
    include("../examples/runexamples.jl")    

end

end