using Test

@testset verbose=true "Modia examples with simulation" begin

include("SimpleFilters.jl")  
include("FilterCircuit.jl")    
include("CauerLowPassFilter.jl")    
include("Rectifier.jl")   
include("MotorControl.jl")    
include("Pendulum.jl")   
include("ServoSystem.jl")

include("Pendulum3D_1.jl")
include("Pendulum3D_2.jl")
include("Pendulum3D_3.jl")
include("BouncingSphere3D_1.jl")
include("BouncingSphere3D_2.jl")
include("Mobile3D.jl")

end
