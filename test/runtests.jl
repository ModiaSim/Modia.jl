module Runtests

using Modia
using Test 

@time begin
    usePlotPackage("SilentNoPlot")
    include("../examples/runexamples.jl")
    usePreviousPlotPackage()
end

end