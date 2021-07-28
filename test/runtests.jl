module Runtests

using Modia
using Modia.ModiaLang
using Test 

@time begin
    usePlotPackage("SilentNoPlot")
    include("$(Modia.ModiaLang.path)/test/runtests.jl")
    include("../examples/runexamples.jl")
    usePreviousPlotPackage()
end

end