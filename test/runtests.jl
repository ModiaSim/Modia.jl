module Runtests

using Modia
using Test 

@time begin
    usePlotPackage("SilentNoPlot")
    include("$(Modia.ModiaLang.path)/test/runtests.jl")
    include("../examples/runexamples.jl")
    usePreviousPlotPackage()
end

end