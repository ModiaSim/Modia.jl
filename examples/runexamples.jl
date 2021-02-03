using Test

@testset "TinyModia examples with simulation" begin

include("SimpleFilters.jl")  

include("FilterCircuit.jl")    

include("MotorControl.jl")    

include("ServoSystem.jl")

end
