module Runtests

# Run all tests with SilentNoPlot (so not plots)

using TinyModia
using Test 

@time @testset verbose=true "TinyModia (with SilentNoPlot)" begin
    usePlotPackage("SilentNoPlot")
    include("include_all.jl")  
    usePreviousPlotPackage()
end

end