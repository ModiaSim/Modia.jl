module Runtests

# Run all tests with SilentNoPlot (so not plots)

using ModiaLang
using Test 

const  test_title = "Test ModiaLang (version=$(ModiaLang.Version) with SilentNoPlot)"

println("\n... $test_title")

@time @testset verbose=true "$test_title" begin
    usePlotPackage("SilentNoPlot")
    include("include_all.jl")  
    usePreviousPlotPackage()
end

end