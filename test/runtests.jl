module Runtests

# Run all tests with SilentNoPlot (so not plots)

using ModiaLang
using Test 

@time @testset verbose=true "ModiaLang (with SilentNoPlot)   " begin
    usePlotPackage("SilentNoPlot")
    include("include_all.jl")  
    usePreviousPlotPackage()
end

end