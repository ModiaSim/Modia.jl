using Test

@testset verbose=true "Modia examples with simulation" begin

include("SimpleFilters.jl")  
include("FilterCircuit.jl")    
include("CauerLowPassFilter.jl")    
include("Rectifier.jl")   
@test_skip include("MotorControl.jl")     # skipped due to issue https://github.com/SciML/DifferentialEquations.jl/issues/950
@test_skip include("Pendulum.jl")         # skipped due to issue https://github.com/SciML/DifferentialEquations.jl/issues/950
include("ServoSystem.jl")

end
