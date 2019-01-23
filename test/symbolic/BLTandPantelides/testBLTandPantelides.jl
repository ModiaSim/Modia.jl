"""
Module with tests for BLTandPantelides.

* Author: Hilding Elmqvist, Mogram AB  
* Date: July-August 2016. Array version: Summer 2018
* License: MIT

"""
module testBLTandPantelides

using Modia.BLTandPantelides
using Modia.BLTandPantelidesUtilities
@static if VERSION < v"0.7.0-DEV.2005"
  using Base.Test
else
  using Test
end

@testset "BLTandPantelides" begin

# testMatch
  println("\nTest match")
  G = Any[
    [2],
    [3,8],
    [4,7],
    [5],
    [3,6],
    [],
     [4,6],
    [1,7],
  ]

  assign = matching(G, 8)
  
  @show assign
  @testset "testMatch" begin
    @test any(assign .== 6) == false
    @test assign == [8,1,2,7,4,5,3,0]
  end


# testSingular
  println("\nSingular system")
  G = Any[
    [3],
    [3],
    [2,3]
  ]

  assign = matching(G, 4)
  @show assign
  
  (invAssign, unAssignedVariables) = invertAssign(assign, length(G))
  @show invAssign, unAssignedVariables

  (ass, unAssignedEquations) = invertAssign(invAssign, length(assign))
  @show ass, unAssignedEquations
  
    @testset "Singular system" begin
    @test assign == [0,3,1,0]
    @test invAssign == [3,0,2]
    @test unAssignedVariables == [1,4]
    @test ass == assign
    @test unAssignedEquations == [2]
  end

  

# testTarjan
  println("\nTest Tarjans strong connect")
  # Reference: Tarjan, Fig.3., page 158.
  G = Any[
    [2],
    [3,8],
    [4,7],
    [5],
    [3,6],
    [],
     [4,6],
    [1,7],
  ]

  assign = 1:8
  
  components = BLT(G, assign)
  @show components
  @testset "testTarjan" begin
    @test components == Any[Any[6],Any[7,5,4,3],Any[8,2,1]]
  end

  

# testPendulum
  println("\nFixed-length pendulum")
  # Reference: Pantelides, page 222-225.
  G = Any[
    [3, 5],
    [4, 6],
    [1, 7, 9],
    [2, 8, 9],
    [1, 2]
  ]

  assign = matching(G, 9)
  @show assign
    
  # Variable vector is: [x, y, w, z, der(x), der(y), der(w), der(z), T], i.e. vector A is:
  A = [5,6,7,8,0,0,0,0,0]

  # Symbolic output of results  

  variables = ["x", "y", "w", "z", "der(x)", "der(y)", "der(w)", "der(z)", "T"]
  equations = [
    "der(x) = w", 
    "der(y) = z",
    "der(w) = T*x",
    "der(z) = T*y - g",
    "0 = x^2 + y^2 - L^2"]
  
  println("\nAssigned original equations:")
  printAssignedEquations(equations, variables, 1:length(equations), assign, A, fill(0, length(equations)))  
  printUnassigned(equations, variables, assign, A, fill(0, length(equations)))
  
  println("\nTest diagnostics for too many equations")
  tooManyEquations = [equations; [
    "x = cos(phi)",
    "y = sin(phi)"]]
  tooFewVariables = [variables; ["phi", "dummy"]]
#=
  Gbig = {G; {
    [1, 10]
    [2, 10]
    }
  }
=#  
  Gbig = copy(G)
  push!(Gbig, [1, 10])
  push!(Gbig, [2, 10])
  @show Gbig  
  Abig = [A; fill(0, 1) ]
  EGbig = [Gbig; buildExtendedSystem(Abig)]
  EGbig = addDependencies(EGbig, (length(Abig)+1):length(EGbig))
  @show EGbig
  equationsBig = [tooManyEquations; fill("h(., der(.)) = 0", length(EGbig)-length(tooManyEquations))]
  assignBig = matching(EGbig, length(EGbig))
  Abig = [Abig; fill(0, length(EGbig)-length(Abig))]
  printAssignedEquations(equationsBig, tooFewVariables, 1:length(EGbig), assignBig, Abig, fill(0, length(EGbig)))  
  printUnassigned(equationsBig, tooFewVariables, assignBig, Abig, fill(0, length(EGbig)))
  componentsBig = BLT(EGbig, assignBig)
  @show componentsBig
  
  println("\nTest diagnostics for too many variables")
  tooManyVariables = [variables; ["dummy"]]
  Gbig = copy(G)
  Gbig[5] = [1, 10] # Wrong equation involving x and dummy
  Gbig[4] = [10, 8, 9] # Wrong equation involving dummy
  @show Gbig  
  Abig = [A; fill(0, 1) ]
  EGbig = [Gbig; buildExtendedSystem(Abig)]
  EGbig = [EGbig; buildFullIncidence(length(Abig)-length(EGbig), length(Abig))]
  @show EGbig
  equationsBig = copy(equations)
  equationsBig[5] = "0 = x^2 + dummy^2 - L^2"
  equationsBig = [equationsBig; fill("h(., der(.)) = 0", length(EGbig)-length(equationsBig))]
  assignBig = matching(EGbig, length(EGbig))
  Abig = [Abig; fill(0, length(EGbig)-length(Abig))]
  printAssignedEquations(equationsBig, tooManyVariables, 1:length(EGbig), assignBig, Abig, fill(0, length(EGbig)))  
  printUnassigned(equationsBig, tooManyVariables, assignBig, Abig, fill(0, length(EGbig)))
  componentsBig = BLT(EGbig, assignBig)
  @show componentsBig

  println("\nTest diagnostics for too few equations")
  tooManyVariables = variables
  Gbig = copy(G)
  pop!(Gbig)  # Delete last equation
  @show Gbig  
  Abig = A
  EGbig = [Gbig; buildExtendedSystem(Abig)]
  equationsBig = equations[1:4]
  equationsBig = [equationsBig; fill("h(., der(.)) = 0", length(EGbig)-length(equationsBig))]
  EGbig = [EGbig; buildFullIncidence(length(Abig)-length(EGbig), length(Abig))]
  equationsBig = [equationsBig; fill("full(...) = 0", length(EGbig)-length(equationsBig))]
  @show EGbig
  assignBig = matching(EGbig, length(EGbig))
  Abig = [Abig; fill(0, length(EGbig)-length(Abig))]
  printAssignedEquations(equationsBig, tooManyVariables, 1:length(EGbig), assignBig, Abig, fill(0, length(EGbig)))  
  printUnassigned(equationsBig, tooManyVariables, assignBig, Abig, fill(0, length(EGbig)))
  componentsBig = BLT(EGbig, assignBig)
  @show componentsBig

  

  println("\nCheck consistency of equations by matching extended equation set")
  EG = [G; buildExtendedSystem(A)]
  @show EG
  assign = matching(EG, 9)
  @show assign
  @testset "testPendulum 1" begin
    @test all(assign .> 0)
  end

  println("\nPerform index reduction")
  (assign, A, B) = pantelides!(G, 9, A)
  @show G

  @show assign
  @show A
  @show B
  @testset "testPendulum 2" begin
  # According to Pantelides, page 224-225.
    @test assign == [0,0,0,0,1,2,7,4,3,9,8]
    @test A == [5,6,7,8,10,11,0,0,0,0,0]
    @test B == [7,8,0,0,6,9,0,0,0]
  end
  
  println("------------------------------------------------------")
  println()
  vActive = fill(true, length(A))
  vActive[[1,3]] .= false
  @show vActive
  assign = matching(G, length(A), vActive)
  @show assign

  components = BLT(G, assign)
  @show components
  @testset "testPendulum BLT components" begin
    @test assign == [0, 5, 0, 2, 1, 6, 7, 4, 3, 9, 8]
    @test components == Any[Any[1], Any[5], Any[6], Any[2], Any[4, 8, 9, 7, 3]]
  end
  
  println("------------------------------------------------------")
  println()
    
  println("\nAll unknowns:")
  printList(variables, 1:length(A), A)
  println("\nAll equations:")
  printList(equations, 1:length(B), B, true)
  println("\nAssigned equations:")
  printAssignedEquations(equations, variables, 1:length(B), assign, A, B)
  println("\nSorted equations:")
  printSortedEquations(equations, variables, components, assign, A, B)
  printUnassigned(equations, variables, assign, A, B)
  
  println("\nBuild augmented system.")
  AG = [G; buildFullIncidence(length(A)-length(B), length(A))]
  @show AG
  assignAG = matching(AG, length(A))
  @show assignAG
  componentsAG = BLT(AG, assignAG)
  @show componentsAG  
  @testset "testPendulum 3" begin
  # See Pantelides, page 222, paragraph 1.
    @test componentsAG == Any[Any[11,3,7,9,8,2,10,4,5,6,1]]
  end

  equationsAG = [equations; fill("", length(B)-length(equations)); fill("full", length(A)-length(B))]
  BAG = [B;fill(0, length(A)-length(B))]
  println("\nAssigned augmented equations:")
  printAssignedEquations(equationsAG, variables, 1:length(BAG), assignAG, A, BAG)
  println("\nSorted augmented equations:")
  printSortedEquations(equationsAG, variables, componentsAG, assignAG, A, BAG)
  printUnassigned(equationsAG, variables, assignAG, A, BAG)
  
  println("\nSet initial conditions on x and y. Should fail.")
  IG1 = copy(G)
  push!(IG1, [1])
  push!(IG1, [2])
  @show IG1
  assignIG1 = matching(IG1, length(A))
  @show assignIG1
  @testset "testPendulum 4" begin
    @test any(assignIG1 .== 0)
  end
  componentsIG1 = BLT(IG1, assignIG1)
  @show componentsIG1

  equationsIG = [equations; fill("", length(B)-length(equations)); fill("initial", length(A)-length(B))]
  BIG = [B;fill(0, length(A)-length(B))]

  printUnassigned(equationsIG, variables, assignIG1, A, BIG)

  println("\nSet initial conditions on x and w.")
  IG2 = copy(G)
  push!(IG2, [1])
  push!(IG2, [3])
  @show IG2
  assignIG2 = matching(IG2, length(A))
  @show assignIG2
  @testset "testPendulum 5" begin
    @test all(assignIG2 .> 0)
  end  
  componentsIG2 = BLT(IG2, assignIG2)
  @show componentsIG2
  
  println("\nSorted IG2 equations:")
  printSortedEquations(equationsIG, variables, componentsIG2, assignIG2, A, BIG)

  println("\nSet initial conditions on w and z.")
  IG3 = copy(G)
  push!(IG3, [3])
  push!(IG3, [4])
  @show IG3
  assignIG3 = matching(IG3, length(A))
  @show assignIG3
  @testset "testPendulum 6" begin
    @test all(assignIG3 .> 0)
  end  
  componentsIG3 = BLT(IG3, assignIG3)
  @show componentsIG3
  
  println("\nSorted IG3 equations:")
  printSortedEquations(equationsIG, variables, componentsIG3, assignIG3, A, BIG)


# testPendulum1
  println("\nFixed-length pendulum")
  # Reference: Pantelides, page 222-225.
  G = Any[
    [3, 5],
    [4, 6],
    [1, 7, 9],
    [2, 8, 9],
    [1, 2]
  ]
    
  # Variable vector is: [x, y, w, z, der(x), der(y), der(w), der(z), T], i.e. vector A is:
  A = [5,6,7,8,0,0,0,0,0]

  println("\nPerform index reduction")
  (assign, A, B) = pantelides!(G, 9, A)
  @testset "testPendulum 2" begin
  # According to Pantelides, page 224-225.
    @test assign == [0,0,0,0,1,2,7,4,3,9,8]
    @test A == [5,6,7,8,10,11,0,0,0,0,0]
    @test B == [7,8,0,0,6,9,0,0,0]
  end
    
  println("\nSet initial conditions on x and w.")
  IG = copy(G)
  push!(IG, [1])
  push!(IG, [3])
  @show IG
  assignIG = matching(IG, length(A))
  @show assignIG
  componentsIG = BLT(IG, assignIG)
  @show componentsIG
  @testset "testPendulum" begin
    @test all(assignIG .> 0)
    @test componentsIG == Any[Any[11], Any[1], Any[10], Any[5], Any[6], Any[2], Any[7,9,8,4,3]]
  end  


# testReactor
  # Reference: Pantelides, page 225-228.
  println("\nExothermic Reactor Model")
  G = Any[
    [1, 3, 5],
    [2, 4, 5, 6],
    [1, 2, 5],
    [1]
  ]
  
  matching(G, 6)

  A = [3,4,0,0,0,0]

  (assign, A, B) = pantelides!(G, 6, A)
  @show assign
  @show A
  @show B
  components = BLT(G, assign)
  @show components
    
  @testset "testReactor" begin
    @test assign == [0,0,1,7,3,2,8,6]
    @test components == Any[Any[3],Any[1],Any[8],Any[6],Any[7],Any[2],Any[4],Any[5]]    
  end

  println("\n\n----------------------\n")
end
  
println("\n----------------------\n")

function bigTest(G)
  n = length(G)
  if n <= 100
    @show G
  end
  assign = matching(G, n)
  if n <= 100
    @show assign
  end
  components = BLT(G, assign)
  if n <= 100
    @show components
  end
end

const n=5000  # Stack overflow for band and n=10000
const nFull=1000
  
function test()
  println("\nBig tests, n = ", n)
  println("\nBig test: diagonal")
  @static if VERSION < v"0.7.0-DEV.2005"
    G1 = Array{Any}(n)
  else
    G1 = Array{Any}(undef, n)
  end
  for i in 1:n
    G1[i] = [i]
  end
  @time bigTest(G1)

  println("\nBig test: band")
  @static if VERSION < v"0.7.0-DEV.2005"
    G2 = Array{Any}(n)
  else
    G2 = Array{Any}(undef, n)
  end
  for i in 1:n
    G2[i] = [i-1 < 1 ? n : i-1, i, mod(i,n)+1]
  end
  @time bigTest(G2)

  println("\nBig test: full, n=", nFull)
  
   @static if VERSION < v"0.7.0-DEV.2005"
    G3 = Array{Any}(nFull)
  else
    G3 = Array{Any}(undef, nFull)
  end
  for i in 1:nFull
    G3[i] = [i for i in 1:nFull]
  end
  @time bigTest(G3)

end

test()

end

