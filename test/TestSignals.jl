module TestSignals

using Modia
using Modia: ValuesID
using Modia.Test

const log = false

if log
    println("\nIndependent Variable")
end
t1 = [0.0, 1.0, 2.0, 3.0]
t2 = [3.0, 4.0, 5.0, 6.0, 7.0]
t3 = [7.0, 8.0, 9.0]
t  = [t1, t2, t3]
sig = Modia.signal(t, t, ValuesID[ValuesID(1, ())], true, log=log, name="t")
@test typeof(sig) == Vector{Float64}
@test sig == [0.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 6.0, 7.0, 7.0, 8.0, 9.0]


if log
    println("\nInvariant tuple variable: Scalar")
end
ts1 = [(false,  1.0, true),
       (false,  2.0, true),
       (false,  3.0, true),
       (false,  4.0, true)]
ts2 = [(false,  5.0, true),
       (false,  6.0, true),
       (false,  7.0, true),
       (false,  8.0, true),
       (false,  9.0, true)]         
ts3 = [(false, 10.0, true),
       (false, 11.0, true),
       (false, 12.0, true)]        
ts = [ts1, ts2, ts3]
sig = Modia.signal(t, ts, ValuesID[ValuesID(2, missing)], false, log=log, name="ts")
@test typeof(sig) == Vector{Float64}
@test sig == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]


if log
    println("\nSegmented tuple variable: Scalar")
end
ts_segmented = [ts1, ts3]  
sig = Modia.signal(t, ts_segmented, ValuesID[ValuesID(1, 2, missing),
                                    ValuesID(3, 2, missing)], false, log=log, name="ts_segmented")
# sig = [1.0, 2.0, 3.0, 4.0, missing, missing, missing, missing, missing, 10.0, 11.0, 12.0]
@test typeof(sig) == Vector{Union{Missing, Float64}}
@test ismissing.(sig) == Bool[0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0]


if log
    println("\nInvariant tuple variable: Vector")
end
tv1 = [(false, [ 1.1,  1.2,  1.3], true),
       (false, [ 2.1,  2.2,  2.3], true),
       (false, [ 3.1,  3.2,  3.3], true),
       (false, [ 4.1,  4.2,  4.3], true)]
tv2 = [(false, [ 5.1,  5.2,  5.3], true),
       (false, [ 6.1,  6.2,  6.3], true),
       (false, [ 7.1,  7.2,  7.3], true),
       (false, [ 8.1,  8.2,  8.3], true),
       (false, [ 9.1,  9.2,  9.3], true)]         
tv3 = [(false, [10.1, 10.2, 10.3], true),
       (false, [11.1, 11.2, 11.3], true),
       (false, [12.1, 12.2, 12.3], true)]        
tv = [tv1, tv2, tv3]
sig = Modia.signal(t, tv, ValuesID[ValuesID(2, missing)], false, log=log, name="tv")
@test typeof(sig) == Matrix{Float64}
@test sig == [1.1 1.2 1.3; 2.1 2.2 2.3; 3.1 3.2 3.3; 4.1 4.2 4.3; 5.1 5.2 5.3; 6.1 6.2 6.3; 7.1 7.2 7.3; 8.1 8.2 8.3; 9.1 9.2 9.3; 10.1 10.2 10.3; 11.1 11.2 11.3; 12.1 12.2 12.3]


if log
    println("\nSegmented tuple variable: Vector")
end
tv_segmented = [tv1, tv3]  
sig = Modia.signal(t, tv_segmented, ValuesID[ValuesID(1, 2, missing),
                                       ValuesID(3, 2, missing)], false, log=log, name="tv_segmented")
# sig = [1.1 1.2 1.3; 2.1 2.2 2.3; 3.1 3.2 3.3; 4.1 4.2 4.3; missing missing missing; missing missing missing; missing missing missing; missing missing missing; missing missing missing; 10.0 missing missing; 11.0 missing missing; 12.0 missing missing]                                        
@test typeof(sig) == Matrix{Union{Missing, Float64}}
@test ismissing.(sig) == Bool[0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 0 0 0; 0 0 0; 0 0 0]


if log
    println("\nSegmented tuple variable: Varying Size Vectors")
end
tsv1 = [(false, [ 1.1,  1.2,  1.3], true),
        (false, [ 2.1,  2.2,  2.3], true),
        (false, [ 3.1,  3.2,  3.3], true),
        (false, [ 4.1,  4.2,  4.3], true)]       
tsv3 = [(false, [10.0], true),
        (false, [11.0], true),
        (false, [12.0], true)]         
tsv = [tsv1, tsv3]
sig = Modia.signal(t, tsv, ValuesID[ValuesID(1, 2, missing),
                     ValuesID(3, 2, missing)], false, log=log, name="tsv")
# sig = [1.1 1.2 1.3; 2.1 2.2 2.3; 3.1 3.2 3.3; 4.1 4.2 4.3; missing missing missing; missing missing missing; missing missing missing; missing missing missing; missing missing missing; 10.0 missing missing; 11.0 missing missing; 12.0 missing missing]             
@test typeof(sig) == Matrix{Union{Missing, Float64}}
@test ismissing.(sig) == Bool[0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 1]


if log                                               
    println("\nInvariant vector variable: Scalar")
end    
v1 = [[-1.0,  1.11,  1.21,  1.31, 1.12,  1.22,  1.32, -2.0],
      [-1.0,  2.11,  2.21,  2.31, 1.12,  1.22,  1.32, -2.0],
      [-1.0,  3.11,  3.21,  3.31, 1.12,  1.22,  1.32, -2.0],
      [-1.0,  4.11,  4.21,  4.31, 1.12,  1.22,  1.32, -2.0]]
v2 = [[-1.0,  5.11,  5.21,  5.31, 1.12,  1.22,  1.32, -2.0],
      [-1.0,  6.11,  6.21,  6.31, 1.12,  1.22,  1.32, -2.0],
      [-1.0,  7.11,  7.21,  7.31, 1.12,  1.22,  1.32, -2.0],
      [-1.0,  8.11,  8.21,  8.31, 1.12,  1.22,  1.32, -2.0],
      [-1.0,  9.11,  9.21,  9.31, 1.12,  1.22,  1.32, -2.0]]
v3 = [[-1.0, 10.11, 10.21, 10.31, 1.12,  1.22,  1.32, -2.0],
      [-1.0, 11.11, 11.21, 11.31, 1.12,  1.22,  1.32, -2.0],
      [-1.0, 12.11, 12.21, 12.31, 1.12,  1.22,  1.32, -2.0]]   
v = [v1, v2, v3]
sig = Modia.signal(t, v, ValuesID[ValuesID(2, ())], true, log=log, name="v")
@test typeof(sig) == Vector{Float64}
@test sig == [1.11, 2.11, 3.11, 4.11, 5.11, 6.11, 7.11, 8.11, 9.11, 10.11, 11.11, 12.11]


if log
    println("\nInvariant vector variable: Vector")
end
sig = Modia.signal(t, v, ValuesID[ValuesID(2, (3,))], true, log=log, name="v")
@test typeof(sig) == Matrix{Float64}
@test sig == [1.11 1.21 1.31; 2.11 2.21 2.31; 3.11 3.21 3.31; 4.11 4.21 4.31; 5.11 5.21 5.31; 6.11 6.21 6.31; 7.11 7.21 7.31; 8.11 8.21 8.31; 9.11 9.21 9.31; 10.11 10.21 10.31; 11.11 11.21 11.31; 12.11 12.21 12.31]


if log
    println("\nInvariant vector variable: Matrix")
end
sig = Modia.signal(t, v, ValuesID[ValuesID(2, (3,2))], true, log=log, name="v")
@test typeof(sig) == Array{Float64, 3}
@test sig == [1.11 1.21 1.31; 2.11 2.21 2.31; 3.11 3.21 3.31; 4.11 4.21 4.31; 5.11 5.21 5.31; 6.11 6.21 6.31; 7.11 7.21 7.31; 8.11 8.21 8.31; 9.11 9.21 9.31; 10.11 10.21 10.31; 11.11 11.21 11.31; 12.11 12.21 12.31;;; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32]


if log
    println("\nSegmented vector variable: Scalar")
end
v_segmented = [v1, v3]
sig = Modia.signal(t, v_segmented, ValuesID[ValuesID(1, 2, ()),
                             ValuesID(3, 2, ())], true, log=log, name="v_segmented")
# sig = [1.11 1.21 1.31; 2.11 2.21 2.31; 3.11 3.21 3.31; 4.11 4.21 4.31; missing missing missing; missing missing missing; missing missing missing; missing missing missing; missing missing missing; 10.11 10.21 10.31; 11.11 11.21 11.31; 12.11 12.21 12.31]
@test typeof(sig) == Vector{Union{Missing, Float64}}
@test ismissing.(sig) == Bool[0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0]


if log
    println("\nSegmented vector variable: Vector")
end    
sig = Modia.signal(t, v_segmented, ValuesID[ValuesID(1, 2, (3,)),
                             ValuesID(3, 2, (3,))], true, log=log, name="v_segmented")
# sig = [1.11 1.21 1.31; 2.11 2.21 2.31; 3.11 3.21 3.31; 4.11 4.21 4.31; missing missing missing; missing missing missing; missing missing missing; missing missing missing; missing missing missing; 10.11 10.21 10.31; 11.11 11.21 11.31; 12.11 12.21 12.31]                             
@test typeof(sig) == Matrix{Union{Missing, Float64}}
@test ismissing.(sig) == Bool[0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 0 0 0; 0 0 0; 0 0 0]


if log
    println("\nSegmented vector variable: Matrix")
end
sig = Modia.signal(t, v_segmented, ValuesID[ValuesID(1, 2, (3,2)),
                                      ValuesID(3, 2, (3,2))], true, log=log, name="v_segmented") 
# sig = [1.11 1.21 1.31; 2.11 2.21 2.31; 3.11 3.21 3.31; 4.11 4.21 4.31; missing missing missing; missing missing missing; missing missing missing; missing missing missing; missing missing missing; 10.11 10.21 10.31; 11.11 11.21 11.31; 12.11 12.21 12.31;;; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32; missing missing missing; missing missing missing; missing missing missing; missing missing missing; missing missing missing; 1.12 1.22 1.32; 1.12 1.22 1.32; 1.12 1.22 1.32]                                      
@test typeof(sig) == Array{Union{Missing, Float64}, 3}
@test ismissing.(sig) == Bool[0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 0 0 0; 0 0 0; 0 0 0;;; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 0 0 0; 0 0 0; 0 0 0]


if log
    println("\nSegmented vector variable: Varying Size Vectors")
end
vv1 = [[-1.0,  1.11,  1.21,  1.31, 1.12,  1.22,  1.32, -2.0],
       [-1.0,  2.11,  2.21,  2.31, 1.12,  1.22,  1.32, -2.0],
       [-1.0,  3.11,  3.21,  3.31, 1.12,  1.22,  1.32, -2.0],
       [-1.0,  4.11,  4.21,  4.31, 1.12,  1.22,  1.32, -2.0]]
vv3 = [[-1.0, -3.0, 10.11, -2.0],
       [-1.0, -3.0, 11.11, -2.0],
       [-1.0, -3.0, 12.11, -2.0]]
vv = [vv1, vv3]       
sig = Modia.signal(t, vv, ValuesID[ValuesID(1, 2, (3,)),
                             ValuesID(3, 3, (1,))], true, log=log, name="vv") 
# sig = [1.11 1.21 1.31; 2.11 2.21 2.31; 3.11 3.21 3.31; 4.11 4.21 4.31; missing missing missing; missing missing missing; missing missing missing; missing missing missing; missing missing missing; 10.11 missing missing; 11.11 missing missing; 12.11 missing missing]
@test typeof(sig) == Matrix{Union{Missing, Float64}}
@test ismissing.(sig) == Bool[0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 1 1 1; 0 1 1; 0 1 1; 0 1 1]
                             
end