module Runtests_withPlot

using ModiaLang
using Test

const  test_title = "Test ModiaLang (version=$(ModiaLang.Version) with " * currentPlotPackage() * ")"

println("\n... $test_title")

@testset "$test_title" begin
    include("include_all.jl")
end

end