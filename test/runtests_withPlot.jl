module Runtests_withPlot

using ModiaLang
using Test

const  test_title = "Test ModiaLang with simulation (with " * currentPlotPackage() * ")"

@testset "$test_title" begin
    include("include_all.jl")
end

end