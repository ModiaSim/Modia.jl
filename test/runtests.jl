module Runtests

using Test

@testset "Test Modia" begin
    
    include("../examples/CauerLowPassFilter.jl")    
    include("../examples/FilterCircuit.jl")   
    include("../examples/MotorControl.jl")          
    include("../examples/Pendulum.jl")    
    include("../examples/ServoSystem.jl")   
    include("../examples/StateSpace.jl")       
    include("../examples/Tutorial.jl")  
    
end

end