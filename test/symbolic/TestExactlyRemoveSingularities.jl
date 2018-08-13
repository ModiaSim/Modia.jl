"""
    module TestExactlyRemoveSingularities - Test Modia/src/symbolic/ExactlyRemoveSingularities.jl
"""
module TestExactlyRemoveSingularities

@static if VERSION < v"0.7.0-DEV.2005"
  using Base.Test
else
  using Test
  using SparseArrays
end

include("../../src/symbolic/ExactlyRemoveSingularities.jl")

using .ExactlyRemoveSingularities

maxNorm(A::AbstractArray{Int}) = maximum(abs.(A))


@testset "Test ExactlyRemoveSingularities" begin

@testset "... upperTrapezoidal(A1,b1)" begin
#println("\n... Test upperTrapezoidal function:----------------------------")
A1 = [1 2 3 4;
      2 2 3 4;
      3 3 3 4;
      9 8 7 6]
b1 = [10, 11, 13, 30]
x1 = [1, 1, 1, 1]

#printobj("A1", A1)
#printobj("b1", b1)
(A1fac, rk1, p11, p12, b1fac) = upperTrapezoidal(A1,b1)
#printobj("A1fac", A1fac)
#printobj("b1fac", b1fac)
#printobj("p11", p11)
#printobj("p12", p12)
#println("\nrk1 = $rk1")

err = maxNorm(A1fac*x1[p12] - b1fac)
@test err == 0
@test rk1 == 4
end


@testset "... upperTrapezoidal(A2,b2)" begin
AA = [1 2 0 4 5 6 7 8
      0 2 0 4 5 6 7 8
      0 0 0 0 0 0 0 0
      0 0 0 0 0 0 0 0      
      0 0 0 2 5 6 7 8
      0 0 0 4 5 6 7 8]
X  = [2 3 4 5 6 7 8 9]'
BB = AA*X
      
C  = [1 2 3 4 5 6
      2 1 0 0 0 0 
      3 0 1 0 0 0
      4 0 0 1 0 0
      5 0 0 0 1 0
      6 0 0 0 0 1]
A2 = C*AA
B2 = C*BB

#printobj("A2", A2)
#printobj("B2", B2)
(A2fac, rk2, p21, p22, B2fac) = upperTrapezoidal(A2,B2)
#printobj("A2fac", A2fac)
#printobj("B2fac", B2fac)
#printobj("p21", p21)
#printobj("p22", p22)
#println("rk2 = $rk2")
err = maxNorm(A2fac*X[p22,:] - B2fac)
#println("norm(A2fac*X[p22,:] - B2fac) = $err")

@test err == 0
@test rk2 == 4
end


@testset "... upperTrapezoidal(A3,b3)" begin

A3 = [ 5 10  15  20;
      -1 -6 -19 -16;
       1  5  15  19;
       5  6  -1 -12;
       4  9  16  29]
x3 = [ 1  2  3  4]'
b3 = A3*x3

#printobj("A3", A3)
#printobj("b3", b3)
(A3fac, rk3, p31, p32, b3fac) = upperTrapezoidal(A3,b3)
#printobj("A3fac", A3fac)
#printobj("b3fac", b3fac)
#printobj("p31", p31)
#printobj("p32", p32)
#println("\nrk3 = $rk3")
err = maxNorm(A3fac*x3[p32,:] - b3fac)
#println("norm(A3fac*x3[p32,:] - b3fac) = $err")

@test err == 0
@test rk3 == 3
end


@testset "... upperTrapezoidal(A4,b4) (dense and sparse form)" begin

A4  = [    1    -1    -1    0      0       0;
           0     0     0    1      1       0;
           1    -1    -1   -1     -1       0;            
           1    -1     0    0      0      -1;
           0     0     0   -1     -1       0]
B4 = fill(0,5,0)

(A4fac, rk4, p41, p42, b4fac) = upperTrapezoidal(A4,B4)

#printobj("A4", A4)
#printobj("A4fac", A4fac)
#printobj("b4fac", b4fac)
#printobj("p41", p41)
#printobj("p42", p42)
#println("\nrk4 = $rk4")

x4 = [5,2,1,3,-1,2]

err = maxNorm(A4fac*x4)

@test err == 0
@test rk4 == 3


#println("... Test sparse algorithm")
A5 = sparse(A4)
B5 = sparse(B4)
(A5fac, rk5, p51, p52, b5fac) = upperTrapezoidal(A5,B5)

#printobj("A5", A5)
#printobj("A5Full", full(A5))

#printobj("A54fac", A5fac)
#printobj("A54facFull", full(A5fac))

#printobj("b54fac", b5fac)
#printobj("p51", p51)
#printobj("p52", p52)
#println("\nrk5 = $rk5")

err = maxNorm(A5fac*x4)

@test err == 0
@test rk5 == 3
          
end

@testset "... Resistor in parallel to a resistor without ground (dense and sparse form)" begin

#   R1.v = R1.p.v - R1.n.v
#      0 = R1.p.i + R1.n.i
#   R2.v = R2.p.v - R2.n.v
#      0 = R2.p.i + R2.n.i
# R1.p.v = R2.p.v
# R1.n.v = R2.n.v
#      0 = R1.p.i + R2.p.i
#      0 = R1.n.i + R2.n.i
#
# After alias elemination
#   R1.v = R1.p.v - R1.n.v
#      0 = R1.p.i + R1.n.i
#   R2.v = R1.p.v - R1.n.v
#      0 = -R1.p.i - R1.n.i

names = ["R1.p.v"; "R1.n.v"; "R1.v"; "R1.p.i"; "R1.n.i"; "R2.v"]

#     R1.p.v R1.n.v R1.v R1.p.i R1.n.i R2.v
Av = [    1    -1    -1    0      0       0;
          0     0     0    1      1       0;
          1    -1    -1   -1     -1       0;     # redundant equation, introduced additionally
          1    -1     0    0      0      -1;
          0     0     0   -1     -1       0];

  @static if VERSION < v"0.7.0-DEV.2005"
    ix = Array{Int64}(0)
  else
    ix = Array{Int64}(undef, 0)
  end
iy = [1,2]
r  = removeSingularities(Av,ix,iy)
# printRemoveSingularities(Av,r,names)

(iya, eqr, ix1, ix2, eqx, A1, A2) = r

@test iya[1] == 2       # R1.n.v can be arbitrarily set
@test eqr    == [3,5]   # Equations 3 and 5 can be removed since redundant



# println("\nResistors with sparse description")
sparse_Av = sparse(Av)
r = removeSingularities(sparse_Av, ix, iy)
#printRemoveSingularities(sparse_Av,r,names)

(iya, eqr, ix1, ix2, eqx, A1, A2) = r

@test iya[1] == 2       # R1.n.v can be arbitrarily set
@test eqr    == [3,5]   # Equations 3 and 5 can be removed since redundant

end


@testset "... Capacitor in parallel to a capacitor with ground (dense and sparse form)" begin

# println("\n\nCapacitor in parallel to a capacitor with ground")
# After alias elemination
#   C1.v = C1.p.v - C1.n.v
#      0 = C1.p.i + C1.n.i
#   C2.v = C1.p.v - C1.n.v
#      0 = -C1.p.i - C1.n.i + ground.i
#      0 = C1.p.v

names = ["C1.p.v"; "C1.n.v"; "C1.v"; "C1.p.i"; "C1.n.i"; "C2.v"; "ground.i"]

#     C1.p.v C1.n.v C1.v C1.p.i C1.n.i C2.v ground.i
Av = [    1    -1    -1    0      0       0     0;
          0     0     0    1      1       0     0;
          1    -1     0    0      0      -1     0;
          0     0     0   -1     -1       0     1;
          1     0     0    0      0       0     0];

ix = [3,6]
iy = [1,2]
r  = removeSingularities(Av,ix,iy);
# printRemoveSingularities(Av,r,names)

(iya, eqr, ix1, ix2, eqx, A1, A2) = r

@test ix1     == [3]   # Dependent state variables: C1.v
@test ix2     == [6]   # Independent state variables: C2.v
@test eqx     == [3]   # Replace equation 3 by A1*x1 + A2*x2 = 0
@test A1[1,1] == 1     # A1 = [1]
@test A2[1,1] == -1    # A2 = [-1]
      
# println("\nCapacitors with sparse description")
sparse_Av = sparse(Av)
r = removeSingularities(sparse_Av,ix,iy);
#printRemoveSingularities(sparse_Av,r,names)

(iya, eqr, ix1, ix2, eqx, A1, A2) = r

@test ix1     == [3]   # Dependent state variables: C1.v
@test ix2     == [6]   # Independent state variables: C2.v
@test eqx     == [3]   # Replace equation 3 by A1*x1 + A2*x2 = 0
@test A1[1,1] == 1     # A1 = [1]
@test A2[1,1] == -1    # A2 = [-1]

end


@testset "... Fourbar kinematic loop (dense and sparse form)" begin

#println("\n\nFourbar kinematic loop")
# R1 = R41 + R4
# R2 = R12 + R1
# R3 = R23 + R2
# R4 = R34 + R3

names = ["R1"; "R2"; "R3"; "R41"; "R12"; "R23"; "R34"]

  #      R1  R2 R3 R41 R12 R23 R34"
Av = [   -1  0  0   1   0   0   0
          1  -1 0   0   1   0   0
          0  1  -1  0   0   1   0
          0  0  1   0   0   0   1]

ix = [4,5,6,7]
iy = fill(0,0)
r  = removeSingularities(Av,ix,iy);
# printRemoveSingularities(Av,r,names)

(iya, eqr, ix1, ix2, eqx, A1, A2) = r

@test ix1      == [4]        # Dependent state variables: R41
@test ix2      == [5,6,7]    # Independent state variables: R12, R23, R34
@test eqx      == [4]        # Replace equation 4 by A1*x1 + A2*x2 = 0
@test A1[1,1]  == -1         # A1 = [-1]
@test A2       == [-1 -1 -1] # A2 = [-1 -1 -1]

sparse_Av = sparse(Av)
r  = removeSingularities(Av,ix,iy);
# printRemoveSingularities(Av,r,names)

(iya, eqr, ix1, ix2, eqx, A1, A2) = r

@test ix1      == [4]        # Dependent state variables: R41
@test ix2      == [5,6,7]    # Independent state variables: R12, R23, R34
@test eqx      == [4]        # Replace equation 4 by A1*x1 + A2*x2 = 0
@test A1[1,1]  == -1         # A1 = [-1]
@test A2       == [-1 -1 -1] # A2 = [-1 -1 -1]

end
end

end
