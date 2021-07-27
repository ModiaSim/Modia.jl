import Test

Test.@testset "Test basic functionality" begin
    include("TestVariables.jl") 
    include("TestFirstOrder.jl") 
    include("TestFirstOrder2.jl")  
    include("TestSource.jl")         
    include("TestStateSelection.jl") 
    include("TestFilterCircuit.jl")
end

Test.@testset "Test units, uncertainties" begin
    include("TestUnits.jl")  
    include("TestUncertainties.jl")
    include("TestUnitsAndUncertainties.jl")
    
    include("TestTwoInertiasAndIdealGear.jl")   
    include("TestTwoInertiasAndIdealGearWithUnits.jl")   
    include("TestTwoInertiasAndIdealGearWithUnitsAndUncertainties.jl")   
    include("TestTwoInertiasAndIdealGearWithMonteCarlo.jl")   
    Test.@test_skip include("TestTwoInertiasAndIdealGearWithUnitsAndMonteCarlo.jl")  # MonteCarlo and Unitful do not yet work together
end


Test.@testset "Test model components" begin
    include("TestSingularLRRL.jl")  
    include("TestStateSpace.jl")  
    include("TestParameter.jl") 
    include("TestHeatTransfer.jl")  
end

Test.@testset "Test events, etc." begin
    include("TestSimpleStateEvents.jl") 
    include("TestSynchronous.jl") 
    include("TestInputOutput.jl")      
    include("TestPathPlanning.jl")  
    include("TestExtraSimulateKeywordArguments.jl")
    Test.@test_skip include("TestBouncingBall.jl")    
end

Test.@testset "Test examples" begin
    include("../examples/runexamples.jl")    
end
