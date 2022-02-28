module Runtests_withPlot

using Modia
using Modia.Test

const  test_title = "Test Modia (version=$(Modia.Version) with " * currentPlotPackage() * ")"

println("\n... $test_title")

@testset "$test_title" begin
    include("include_all.jl")
end

end