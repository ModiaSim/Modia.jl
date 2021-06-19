module Runtests_withPlot

using TinyModia
using Test

const  test_title = "Test TinyModia with simulation (with " * currentPlotPackage() * ")"

@testset "$test_title" begin
    include("include_all.jl")
end

end