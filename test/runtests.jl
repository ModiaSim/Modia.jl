module Runtests

# Run all tests with SilentNoPlot (so not plots)

using TinyModia
using Test 

@testset "Test TinyModia with simulation (without plots)" begin
    usePlotPackage("SilentNoPlot")
    include("include_all.jl")  
    usePreviousPlotPackage()
end

end