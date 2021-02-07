using Test

@testset "Test TinyModia with simulation" begin

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

include("TestCauerLowPassFilter.jl")    

include("TestSingularLRRL.jl")  

include("../examples/runexamples.jl")    

end
