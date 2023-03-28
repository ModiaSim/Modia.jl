import Modia
import Modia.Test

Test.@testset "Test basic functionality" begin
    include("TestVariables.jl")
    include("TestFirstOrder.jl")
    include("TestFirstOrder2.jl")
    @test_skip include("TestPendulum.jl")    # skipped due to issue https://github.com/SciML/DifferentialEquations.jl/issues/950
    include("TestSource.jl")
    include("TestLinearEquations.jl")
    include("TestStateSelection.jl")
    include("TestFilterCircuit.jl")
    include("TestFilterCircuit2.jl")
    include("TestArrays.jl")
    include("TestLinearSystems.jl")
end


Test.@testset "Test units, uncertainties" begin
    include("TestUnitAsString.jl")
    include("TestUnits.jl")
    @test_skip include("TestUncertainties.jl")              # skipped due to issue https://github.com/SciML/DifferentialEquations.jl/issues/950
    @test_skip include("TestUnitsAndUncertainties.jl")      # skipped due to issue https://github.com/SciML/DifferentialEquations.jl/issues/950
    @test_skip include("TestTwoInertiasAndIdealGear.jl")    # skipped due to issue https://github.com/SciML/DifferentialEquations.jl/issues/950
    include("TestTwoInertiasAndIdealGearWithUnits.jl")
    @test_skip include("TestTwoInertiasAndIdealGearWithUnitsAndUncertainties.jl")    # skipped due to issue https://github.com/SciML/DifferentialEquations.jl/issues/950
    include("TestTwoInertiasAndIdealGearWithMonteCarlo.jl") 
    include("TestTwoInertiasAndIdealGearWithUnitsAndMonteCarlo.jl")
    include("TestLinearEquationSystemWithUnitsAndMonteCarlo.jl")
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
    #Test.@test_skip include("TestBouncingBall.jl")
end


Test.@testset "Test multi returning functions" begin
    include("TestMultiReturningFunction.jl")
    include("TestMultiReturningFunction4A.jl")
    include("TestMultiReturningFunction5A.jl")
    include("TestMultiReturningFunction6.jl")
    include("TestMultiReturningFunction7A.jl")
    include("TestMultiReturningFunction10.jl")
end

include("../examples/runexamples.jl")

