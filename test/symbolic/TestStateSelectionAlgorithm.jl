################################################
#
# Module to test Modia/src/symbolic/StateSelection.jl
#
# Author: Martin Otter, DLR-SR
#
################################################
module TestStateSelectionAlgorithm

@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

include(joinpath("..", "..", "src", "symbolic", "StateSelection.jl"))
using .StateSelection

const printSortedEquations = true

isAlgebraic(v::Int, A::Vector{Int}, Arev::Vector{Int}) = A[v] == 0 && Arev[v] == 0

"""
    checkSortedEquationGraph(eq,Vx,Vderx) - check result of sortedEquationGraph(..)
    
- eq   : Result of sortedEquationGraph(..)
- Vx   : Required value of eq.Vx
- Vderx: Required value of eq.Vderx
"""
function checkSortedEquationGraph(eq;nx=0, nc=0, Vx=Int[], Vderx=Int[])
    @test length(eq.Vx) + length(eq.Vmue) == nx
    @test length(eq.Er)                 == nx
    @test eq.nc    == nc
    @test eq.Vx    == Vx
    @test eq.Vderx == Vderx
   
   # Check correct indices
    for i in eachindex(eq.Vx)
        vx = eq.Vx[i]
        if vx != 0
            if eq.Arev[vx] > 0
                vxInt_index_isZero = eq.VxRev[ eq.Arev[vx] ] == 0
                @test vxInt_index_isZero == false
            end
        end    
    end
   
    for i in eachindex(eq.Vderx)
        vderx = eq.Vderx[i]
        if vderx == 0
            if !isAlgebraic(eq.Vx[i], eq.A, eq.Arev)
                vx = eq.Vx[i]
                dervx = eq.A[vx]
                @test (dervx == 0) == false                            
                if dervx != 0
                    der_vx_index = eq.VxRev[ dervx ]
                    @test (der_vx_index == 0) == false
                end
            end
        end
    end
end



@testset "\nTest StateSelection algorithm" begin 

    @testset "... Test two coupled inertias (all/subset/no unknowns can be solved for)" begin

    #= -----------------------------------------
    @model CoupledInertias1 begin
        J1_J = 1
        J1_tau = Float()
        J1_phi = Float(start = 1.0)
        J1_w = Float(start = 0.0)

        J2_J = 2
        J2_tau = Float()
        J2_phi = Float(start = 0.0)
        J2_w = Float(start = 0.0)
        @equations begin
            J1_w = der(J1_phi)
            J1_J * der(J1_w) = J1_tau

            J2_w = der(J2_phi)
            J2_J * der(J2_w) = J2_tau

            J2_phi = J1_phi
            J2_tau + J1_tau = 0
        end
    end

gives:

Assigned equations after index reduction:
   1:              der(J1_phi):  J1_w = der(J1_phi)
   2:                der(J1_w):  J1_J * der(J1_w) = J1_tau
   3:              der(J2_phi):  J2_w = der(J2_phi)
   4:                   J2_tau:  J2_J * der(J2_w) = J2_tau
   5:                         :  J2_phi = J1_phi
   6:                   J1_tau:  J2_tau + J1_tau = 0
   7:                         :  DER( J2_phi = J1_phi )
   8:             der2(J1_phi):  DER( J1_w = der(J1_phi) )
   9:                der(J2_w):  DER( J2_w = der(J2_phi) )
  10:             der2(J2_phi):  DER2( J2_phi = J1_phi )

Unassigned variables:
J1_phi, J1_w, J2_phi, J2_w

Unassigned equations:
J2_phi = J1_phi
DER(J2_phi = J1_phi)

NEW STATE SELECTION

vNames:
    1: J1_tau
    2: J1_phi                        A[2] = 7
    3: J1_w                          A[3] = 8
    4: J2_tau
    5: J2_phi                        A[5] = 9
    6: J2_w                          A[6] = 10
    7: der(J1_phi)                   A[7] = 11
    8: der(J1_w)
    9: der(J2_phi)                   A[9] = 12
   10: der(J2_w)
   11: der2(J1_phi)
   12: der2(J2_phi)

=#
        Gorigin = [[3,7],[8,1],[6,9],[10,4],[5,2],[4,1],[9,7],[8,11],[10,12],[12,11]]
        BLT = Array{Int64,1}[[1],[3],[8,10,9,4,6,2],[5],[7]]
        assign = [6,0,0,4,0,0,1,2,3,9,8,10]
        Avar = [0,7,8,0,9,10,11,0,12,0,0,0]
        Bequ = [8,0,9,0,7,0,10,0,0,0]
        VNames = ["J1_tau",
          "J1_phi",
          "J1_w",
          "J2_tau",
          "J2_phi",
          "J2_w",
          "der(J1_phi)",
          "der(J1_w)",
          "der(J2_phi)",
          "der(J2_w)",
          "der2(J1_phi)",
          "der2(J2_phi)"]

        # First version where all unknowns can be solved for
        println("\n\n... Test two coupled inertias (all unknowns can be solved for)")
        # Gsolvable = Gorigin
        Gsolvable = StateSelection.newRaggedIntMatrix(length(Gorigin))
        Gsolvable[1:6] = Gorigin[1:6]
        eq = getSortedEquationGraph(Gorigin, Gsolvable, BLT, assign, Avar, Bequ, VNames)
        if printSortedEquations
           printSortedEquationGraph(eq)
        end
        checkSortedEquationGraph(eq; nx=2, nc=0, Vx=[2,7], Vderx=[0,11])


        println("\n\n... Test two coupled inertias (only a subset of unknowns can be solved for)")
        Gsolvable = Any[Any[3,7],Any[1],Any[6,9],Any[4],Any[2,5],Any[1,4],Any[],[],[],[]]
        eq = getSortedEquationGraph(Gorigin, Gsolvable, BLT, assign, Avar, Bequ, VNames)
        if printSortedEquations
           printSortedEquationGraph(eq)
        end
        checkSortedEquationGraph(eq; nx=2, nc=0, Vx=[5,9], Vderx=[0,12])

        println("\n\n... Test two coupled inertias (no unknowns can be solved for)")
        Gsolvable = StateSelection.newRaggedIntMatrix(length(Gorigin))
        eq = getSortedEquationGraph(Gorigin, Gsolvable, BLT, assign, Avar, Bequ, VNames)
        if printSortedEquations
           printSortedEquationGraph(eq)
        end
        checkSortedEquationGraph(eq; nx=9, nc=4, Vx=[2, 5, 7, 9, 6, 3, 0, 0], Vderx=[0, 0, 11, 12, 10, 8, 4, 1])
    end



    @testset "... Test simple sliding mass model" begin


#= 

  1: 0 = r - n*s
  2: 0 = v - der(r)
  3: 0 = m*der(v) - f - m*g - u
  4: 0 = n*f
  5: 0 = u + (c*s + d*der(s) + sf)*n
  9: 0 = der(sf) + sf - su(t)     # additional equation to test solving for der(sf)

  Variables
   1: s
   2: r
   3: v
   4: f
   5: u
  ------------
   6: der(s)
   7: der2(s)
   8: der(r)
   9: der2(r)
  10: der(v)
  11: sf
  12: der(sf)

   Equations
   1: (1)
   2: (2)
   3: (3)
   4: (4)
   5: (5)
   6: der(1)
   7: der2(1)
   8: der(2)
   9: (9)
=#

        println("\n\n... Test simple sliding mass model with Tearing")
        assign = [0,0,0,3,5,0,4,0,7,8,0,9]
        A = [6,8,10,0,0,7,0,9,0,0,12,0]
        B = [6,8,0,0,0,7,0,0,0]
        VNames = ["s", "r", "v", "f", "u", "der(s)", "der2(s)",
          "der(r)", "der2(r)", "der(v)", "sf", "der(sf)"]

        BLT = [[5],
       [9],
       [7,8,3,4]]

        G = [[2,1],
     [3,8],
     [10,4,5],
     [4],
     [5,1,6],
     [8,6],
     [9,7],
     [10,9],
     [12,11]]

    #=
        Gsolvable = [[2],
             [3,8],
             [10,4,5],
             [],
             [5],
             [8],
             [9],
             [10,9],
             [12]]
    =#

        Gsolvable = [[2],
             [3,8],
             [10,4,5],
             [],
             [5],
             [],
             [],
             [],
             [12]]

        eq = getSortedEquationGraph(G, Gsolvable, BLT, assign, A, B, VNames)
        if printSortedEquations
           printSortedEquationGraph(eq)
        end
        checkSortedEquationGraph(eq; nx=3, nc=0, Vx=[1,6,11], Vderx=[0,7,12])

    end



    @testset "... Test Multi-index DAE (with and without tearing)" begin

    #= 
  1: 0 = u1(t) + x1 - x2
  2: 0 = u2(t) + x1 + x2 - x3 + der(x6)
  3: 0 = u3(t) + x1 + der(x3) - x4
  4: 0 = u4(t) + 2*der2(x1) + der2(x2) + der2(x3) + der(x4) + der3(x6)
  5: 0 = u5(t) + 3*der2(x1) + 2*der2(x2) + x5
  6: 0 = u6(t) + 2*x6 + x7
  7: 0 = u7(t) + 3*x6 + 4*x7
  8: 0 = u8(t) + x8 - sin(x8)


  Variables
   1: x1
   2: x2
   3: x3
   4: x4
   5: x5          (5)
   6: x6
   7: x7          (6)
   8: der(x1)
   9: der(x2)
  10: der(x3)
  11: der(x4)     (4)
  12: der(x6)
  13: der2(x1)    der2(1)
  14: der2(x2)    der2(2)
  15: der2(x3)    der(3)
  16: der2(x6)
  17: der3(x6)    der3(7)
  18: x8
  ------------
  19: der(x7)     der(6)
  20: der2(x7)    der2(6)
  21: der3(x7)    der3(6)


   Equations
   1: (1)
   2: (2)
   3: (3)
   4: (4)
   5: (5)
   6: (6)
   7: (7)
   8: (8)
   9: der(1)
  10: der2(1)
  11: der(2)
  12: der2(2)
  13: der(3)
  14: der(6)
  15: der(7)
  16: der2(6)
  17: der2(7)
  18: der3(6)
  19: der3(7)
=#


        assign = [0,0,0,0,5,0,6,0,0,0,4,0,10,12,13,0,19,8,14,16,18]
        A = [8,9,10,11,0,12,19,13,14,15,0,16,0,0,0,17,0,0,20,21,0]
        B = [9,11,13,0,0,14,15,0,10,0,12,0,0,16,17,18,19,0,0]
        VNames = ["x1", "x2", "x3", "x4", "x5", "x6", "x7",
          "der(x1)", "der(x2)", "der(x3)", "der(x4)", "der(x6)",
          "der2(x1)", "der2(x2)", "der2(x3)", "der2(x6)",
          "der3(x6)", "x8", "der(x7)", "der2(x7)", "der3(x7)"]

        BLT = [[8],
       [18,19],
       [10,12,13,4,5]]

        G = [[1,2],
     [1,2,3,12],
     [1,10,4],
     [13,14,15,11,17],
     [13,14,5,18],
     [6,7],
     [6,7],
     [18],
     [8,9],
     [13,14],
     [8,9,10,16],
     [13,14,15,17],
     [8,15,11],
     [12,19],
     [12,19],
     [16,20],
     [16,20],
     [17,21],
     [17,21]]

        println("\n\n... Test Multi-Index DAE without tearing")
        Gsolvable = StateSelection.newRaggedIntMatrix(length(G))
        eq = getSortedEquationGraph(G, Gsolvable, BLT, assign, A, B, VNames)
        if printSortedEquations
           printSortedEquationGraph(eq)
        end
        checkSortedEquationGraph(eq; nx=21, nc=12, Vx=[7, 6, 19, 12, 20, 16, 1, 2, 3, 8, 9, 10, 4, 18, 0], 
                                           Vderx=[0, 0, 0, 0, 21, 17, 0, 0, 0, 13, 14, 15, 11, 0, 5])

                                           
        println("\n\n... Test Multi-Index DAE WITH tearing")
        # Gsolvable    = copy(G)
        # Gsolvable[8] = fill(0, 0)  # "0 = u8(t) + x8 - sin(x8)" cannot be solved for x8
        Gsolvable[1:7] = G[1:7]   
        eq = getSortedEquationGraph(G, Gsolvable, BLT, assign, A, B, VNames)
        if printSortedEquations
           printSortedEquationGraph(eq)
        end
        checkSortedEquationGraph(eq; nx=8, nc=4, Vx=[7, 19, 20, 2, 9, 18], 
                                 Vderx=[0, 0, 21, 0, 14, 0])

    end
end
end
