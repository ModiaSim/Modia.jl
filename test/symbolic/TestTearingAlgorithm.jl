################################################
#
# Module to test Modia/src/symbolic/Tearing.jl
#
# Author: Martin Otter, DLR-SR
#
################################################
module TestTearing

@static if VERSION < v"0.7.0-DEV.2005"
  using Base.Test
else
  using Test
end

include(joinpath("..", "..", "src", "symbolic", "Tearing.jl"))

@testset "\nTest Tearing" begin


    @testset "... Test equation system 1" begin

        G = Any[[1,2,3,4,5], [1,2,3,4,8], [4,1,2], [3,2,1], [2,1,7]]
        es = [2,3,5,1]
        vs = [1,3,7,8]
        td = TraverseDAG(G, 8)

        (eSolved, vSolved, eResidue, vTear) = tearEquations!(td, G, es, vs)

        @test vTear    == [3,8]
        @test eSolved  == [2,5]
        @test vSolved  == [1,7]
        @test eResidue == [3,1]
    end



    @testset "... Test equation system 2" begin

   #=
     eSolved  = [2,5]
     vSolved  = [1,7]
     eResidue = [3,1]
     vTear    = [3,8]
     ignored variables: 2,4,5,6

     known  : 2,3,4,5,6,8
     unknown: 1,7
     [1,2,3,4,8] is solved for 1
     [2,1,7]     is solved for 7
   =#


   #=
      unknowns: x1dd, x2dd, x3dd, x4d, x5
       1: 0 = x1dd - x2dd
      2: 0 = x1dd + x2dd - x3dd + x6ddd
      3: 0 = x1d  + x3dd - x4d
      4: 0 = 2*x1dd + x2dd + x3dd + x4d + x6ddd
      5: 0 = 3*x1dd + 2*x2dd + x5

      Variables
      1: x1dd
      2: x2dd
      3: x3dd
      4: x4d
      5: x5
      6: x1d
      7: x6ddd

      vTear    = [2]
      eSolved  = [1,2,3,5]
      vSolved  = [1,3,4,5]
      eResidue = [4]
   =#

        G = Any[[1,2],
           [1,2,3,7],
           [6,3,4],
           [1,2,3,4,7],
           [1,2,5]]
        es = [1,2,3,4,5]
        vs = [1,2,3,4,5]

        td2 = TraverseDAG(G, 7)
        (eSolved, vSolved, eResidue, vTear) = tearEquations!(td2, G, es, vs)

        @test vTear    == [2]
        @test eSolved  == [1,2,3,5]
        @test vSolved  == [1,3,4,5]
        @test eResidue == [4]
    end



    @testset "... Test equation system 3" begin

        G = Any[ [1,2,4,7],
            [6,3,4],
            [2,5],
            [4,5],
            [1,2,3,7]]
        es = [1,2,3,4,5]
        vs = [1,2,3,4,5]
        td3 = TraverseDAG(G, 7)
        (eSolved, vSolved, eResidue, vTear) = tearEquations!(td3, G, es, vs)

        @test vTear    == [5]
        @test eSolved  == [4,3,1,2]
        @test vSolved  == [4,2,1,3]
        @test eResidue == [5]
    end



    @testset "... Test equation system 4" begin
   #=
      Same as test2, but indices changed

      unknowns: z1, z2, z4, z6, z7
      1: 0 = f1(z2,z6,z7)
       2: 0 = f2(z1,z7)
      3: 0 = f3(z1,z2,z3,z7)
      4: 0 = f4(z6,z3,z4)
      5: 0 = f5(z1,z2,z3,z4,z7)
      6: 0 = f6(z1,z2,z5)
      7: 0 = f7(z4,z5,z6)

      Variables
      1: z1    (x1dd)
      2: z2    (x6ddd)
      3: z3    (x3dd)
      4: z4    (x4d)
      5: z5    (x5)
      6: z6    (x1d)
      7: z7    (x2dd)
   =#

        G = Any[[2,6,7],
           [1,2],
           [1,2,3,7],
           [6,3,4],
           [1,2,3,4,7],
           [1,2,5],
           [4,5,6]]
        es = [2,3,4,5,6]
        vs = [1,3,4,5,7]

        td4 = TraverseDAG(G, 7)
        (eSolved, vSolved, eResidue, vTear) = tearEquations!(td4, G, es, vs)

        @test vTear    == [7]
        @test eSolved  == [2,3,4,6]
        @test vSolved  == [1,3,4,5]
        @test eResidue == [5]
    end



    @testset "... Test equation system 5" begin
   #=
     Check "marked" flag

      unknowns: z1, z2, z3, z4, z5, z6
      1: 0 = f1(z1,z6)
       2: 0 = f2(z2,z1)
      3: 0 = f3(z3,z2)
      4: 0 = f4(z4,z3)
      5: 0 = f5(z5,z4)
      6: 0 = f6(z6,z1)
   =#

        G = Any[[1,6],
           [2,1],
           [3,2],
           [4,3],
           [5,4],
           [6,5]]
        es = [1,2,3,4,5,6]
        vs = [1,2,3,4,5,6]

        td5 = TraverseDAG(G, 6)
        (eSolved, vSolved, eResidue, vTear) = tearEquations!(td5, G, es, vs)

        @test vTear    == [6]
        @test eSolved  == [1,2,3,4,5]
        @test vSolved  == [1,2,3,4,5]
        @test eResidue == [6]
    end



    @testset "... Test equation system 6" begin
   #=
     Check predefined equations

     phi1 = phi2
     w1 = der(phi1)
     w2 = der(phi2)
     der(phi1) = der(phi2)

      unknowns: 1:w1, 2:w2, 3:der(phi1), 4:der(phi2)
      1: 0 = f1(w1,der(phi1))
      2: 0 = f2(w2,der(phi2))
      3: 0 = f3(der(phi1), der(phi2))
   =#

        G = Any[[1,3],
           [2,4],
           [3,4]]
        es = [1,2,3]
        vs = [1,2,3,4]

        td = TraverseDAG(G, 4)
        (eSolved, vSolved, eResidue, vTear) = tearEquations!(td, G, es, vs;eSolvedFixed=[3], vSolvedFixed=[3], vTearFixed=[4])

        @test vTear    == [4]
        @test eSolved  == [3,1,2]
        @test vSolved  == [3,1,2]
        @test eResidue == Int[]
    end



    @testset "... Test equation systems 7a and 7b" begin
        neq = 10
        G = fill(fill(0, 0), neq)
        for i in 1:neq
            G[i] = i == 1 ? [1,neq] : [i,i - 1]
        end
   #G = [[i, (i==1 ? neq : i-1)] for i = 1:neq]
        es = collect(1:neq)
        vs = es
        td5 = TraverseDAG(G, neq)
        (eSolved, vSolved, eResidue, vTear) = tearEquations!(td5, G, es, vs)
        @test vTear == [neq]

        for i in 1:neq
            G[i] = i == neq ? [1,neq] : [neq - i + 1,neq - i]
        end
   #println("G = ")
   #display(G)
   #G = [[i, (i==1 ? neq : i-1)] for i = neq:-1:1]
        es = collect(1:neq)
        vs = es
        td5 = TraverseDAG(G, neq)
        (eSolved, vSolved, eResidue, vTear) = tearEquations!(td5, G, es, vs)
        @test vTear == [1]
    end

end  # facts

end  # module
