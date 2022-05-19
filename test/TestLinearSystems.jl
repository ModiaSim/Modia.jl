module TestLinearSystems

using Modia
using Modia.OrderedCollections
using Modia.LinearAlgebra
using Modia.StaticArrays
using Modia.Test
@usingModiaPlot

"""
    ls = LinearStateSpace(; A, B, C, x_init=nothing)

Define a linear state space system Model:

```math
\begin{aligned}
\frac{dx}{dt} &= A*x + B*u \\
            y &= C*x
\end{aligned}
```

where

- `x_init` is the optional vector of init-values. If not provided, zeros will be used.
- The number of **inputs** (= `size(B,2)`) and **outputs** (= `size(C,1)`) is **fixed**, after @instantiateModel(..) was called.
- The number of **states** (= `size(A,1)`) can be **changed** before simulation starts,
  by providing appropriate ``A,B,C`` matrices and `x_init` vector as `merge` values in `simulate!(..., merge = ...)`.

# Example
```
using Modia
@usingModiaPlot

# T*der(x) + x = u
T = 0.2;
SSTest = Model(
            ss = LinearStateSpace(A=[-1.0/T;;], B=[1.0/T;;], C=[0.9;;], x_init=[0.2]), # one state
            equations = :[ss.u = 2.0,
                          y = ss.y[1]]
         )

ssTest = @instantiateModel(SSTest, logCode=false)
simulate!(ssTest, stopTime=1.0, log=true, logStates=true)
plot(ssTest, ("ss.x", "ss.u", "y"), figure=1)

simulate!(ssTest, stopTime=1.0, log=true, logStates=true,
                                merge=Map(ss = Map(A=[-1/T   0.0;
                                                       0.0  -1/T],
                                                   B=[1.0/T;
                                                      1.0/T;;],
                                                   C=[0.4 0.4;],
                                                   x_init=[0.2,0.4]))) # two states
plot(ssTest, ("ss.x", "ss.u", "y"), figure=2)
```
"""
LinearStateSpace(; kwargs...) = Model(; _buildFunction = :(buildLinearStateSpace!),         # Called once in @instantiateModel(..) before getDerivatives!(..) is generated
                                        _instantiateFunction = Par(functionName = :(instantiateLinearStateSpace!)),  # Called once after new A,B,C values are merged
                                        kwargs...)

mutable struct LinearStateSpaceStruct{FloatType}
    path::String              # Path name of instance
    x_hidden_startIndex::Int
    A::Matrix{FloatType}
    B::Matrix{FloatType}
    C::Matrix{FloatType}
    x_init::Vector{FloatType}  # Initial values of states
    y::Vector{FloatType}       # Internal memory for y
    x::Vector{FloatType}       # Internal memory for x
    der_x::Vector{FloatType}   # Internal memory for der(x)

    function LinearStateSpaceStruct{FloatType}(; A::AbstractMatrix, B::AbstractMatrix, C::AbstractMatrix,
                                                 x_init::Union{AbstractVector,Nothing}=nothing,
                                                 u::AbstractVector, y::AbstractVector,  # Code generated with buildLinearStateSpace! provides start values of u and y.
                                                 path::String, kwargs...) where {FloatType}
        #println("... 4: LinearStateSpaceStruct called for $path")
        if length(kwargs) > 0
            @warn "LinearStateSpaceStruct with path=$path ignores keyword arguments: $(kwargs...)"
        end
        @assert(size(A,2) == size(A,1))
        @assert(size(B,1) == size(A,1))
        @assert(size(C,2) == size(A,1))
        @assert(size(u,1) == size(B,2))
        @assert(size(y,1) == size(C,1))
        @assert(size(u,1) == size(B,2))
        @assert(size(y,1) == size(C,1))
        copyA = Matrix{FloatType}(deepcopy(A))
        copyB = Matrix{FloatType}(deepcopy(B))
        copyC = Matrix{FloatType}(deepcopy(C))
        copy_x_init = if isnothing(x_init); zeros(FloatType, size(A,1)) else Vector{FloatType}(deepcopy(x_init)) end
        new(path, 0, copyA, copyB, copyC, copy_x_init, zeros(FloatType,size(C,1)), zeros(FloatType, size(A,1)), zeros(FloatType,size(A,1)))
    end
end

mutable struct LinearStateSpaceBuild{FloatType}
    # Cannot be changed after @instantiatedModel(..)
    path::String # Path name of instance, e.g. "a.b.c"
    nu::Int      # Number of inputs
    ny::Int      # Number of outputs

    # Can be changed after @instantiateModel
    ls::Union{LinearStateSpaceStruct{FloatType}, Nothing}

    function LinearStateSpaceBuild{FloatType}(path::String, nu::Int, ny::Int) where {FloatType}
        #println("... 2: LinearStateSpaceBuild called with path = $path")
        @assert(nu >= 0)
        @assert(ny >= 0)
        new(path,nu,ny,nothing)
    end
end


function buildLinearStateSpace!(model::AbstractDict, FloatType::Type, TimeType::Type,
                                buildDict::OrderedCollections.OrderedDict{String,Any},
                                path::Union{Expr,Symbol,Nothing})
    # Called from @instantiatedModel, during instantiation of the model.
    pathAsString = Modia.modelPathAsString(path)
    #println("... 1: buildLinearStateSpace! called for path = ", pathAsString)

    # Determine nu,ny from model
    B = model[:B]
    C = model[:C]
    nu = size(B,2)
    ny = size(C,1)
    u_zeros = zeros(FloatType,nu)
    y_zeros = zeros(FloatType,ny)

    # Define code to be generated
    lsCode = Model(ls      = Var(hideResult=true),
                   success = Var(hideResult=true),
                         u = Var(input  = true, start = u_zeros),
                         y = Var(output = true, start = y_zeros),
                    equations = :[
                        ls = openLinearStateSpace!(instantiatedModel, $pathAsString)
                        y = computeOutputs!(instantiatedModel, ls)
                        success = computeStateDerivatives!(instantiatedModel, ls, u)])

    # Store build info in buildDict
    buildDict[pathAsString] = LinearStateSpaceBuild{FloatType}(pathAsString, nu, ny)
    return lsCode
end


function instantiateLinearStateSpace!(partiallyInstantiatedModel::SimulationModel{FloatType,TimeType},
                                      model::AbstractDict, path::String; log=false)::Nothing where {FloatType,TimeType}
    # Called during evaluation of the parameters (before initialization)
    if log
        println("instantiateLinearStateSpace! called for $path with model = $model")
    end
    lsBuild::LinearStateSpaceBuild{FloatType} = partiallyInstantiatedModel.buildDict[path]
    ls = LinearStateSpaceStruct{FloatType}(; path, model...)
    @assert(size(ls.A,2) == size(ls.A,1))
    @assert(size(ls.B,2) == lsBuild.nu)
    @assert(size(ls.C,1) == lsBuild.ny)
    ls.x_hidden_startIndex = Modia.addHiddenState!(partiallyInstantiatedModel, path*".x", path*".der(x)", ls.x_init)
    lsBuild.ls = ls
    return nothing
end


function openLinearStateSpace!(instantiatedModel::SimulationModel{FloatType,TimeType}, path::String)::LinearStateSpaceStruct{FloatType} where {FloatType,TimeType}
    ls = instantiatedModel.buildDict[path].ls
    Modia.copyFromHiddenState!(instantiatedModel, ls.x_hidden_startIndex, ls.x)
    return ls
end

function computeOutputs!(instantiatedModel, ls)
    mul!(ls.y, ls.C, ls.x)
    return ls.y
end

function computeStateDerivatives!(instantiatedModel, ls, u)::Bool
    # ls.der_x .= ls.A*ls.x + ls.B*u
    mul!(ls.der_x, ls.A, ls.x)
    mul!(ls.der_x, ls.B, u, 1.0, 1.0)
    Modia.copyToHiddenStateDerivative!(instantiatedModel, ls.x_hidden_startIndex, ls.der_x)
    return true
end

# T*der(x) + x = u
T = 0.2;
SSTest = Model(
            ss = LinearStateSpace(A=[-1.0/T;;], B=[1.0/T;;], C=[0.9;;], x_init=[0.2]),  # one state
            equations = :[ss.u = [2.0],
                          y = ss.y[1]]
         )

ssTest = @instantiateModel(SSTest, logCode=true)
simulate!(ssTest, stopTime=1.0, log=false, logStates=true, requiredFinalStates = [1.987867388853733])
Modia.printResultInfo(ssTest)
plot(ssTest, ("ss.x", "ss.u", "y"), figure=1)

simulate!(ssTest, stopTime=1.0, log=false, logStates=true,
                                merge=Map(ss = Map(A=[-1/T   0.0;
                                                       0.0  -1/T],
                                                   B=[1.0/T;
                                                      1.0/T;;],
                                                   C=[0.4 0.4;],
                                                   x_init=[0.2,0.4])),  # two states
          requiredFinalStates = [1.98786636233743, 1.9892145443000466])
Modia.printResultInfo(ssTest)
plot(ssTest, ("ss.x", "ss.u", "y"), figure=2)


println("\n... Test init vectors of scalars, fixed-size, variable-size, hidden-size vectors")

T = 0.2
                
SSTest2 = Model(
                submodel = Model(   ss = LinearStateSpace(A=[-1.0/T;;], B=[1.0/T;;], C=[0.9;;], x_init=[0.2]), 
                                    x1 = Var(init = 1.1), 
                                    x3 = Var(init = SVector{3}(0.5, 0.6, 0.7)),
                                    x4 = Var(init = [0.8, 0.9]),
                                equations = :[ der(x1) = -x1
                                               der(x2) = -x2
                                               der(x3) = -x3
                                               der(x4) = -x4 ]),
                ss = LinearStateSpace(A=[-1/T   0.0;
                                         0.0   -1/T],
                                      B=[1.0/T;
                                         1.0/T;;],
                                      C=[0.4 0.4;],
                                      x_init=[0.3,0.4]),
            equations = :[ss.u = [2.0], submodel.ss.u = [2.1],
                          y1 = ss.y[1], y2 = submodel.ss.y[1]]                                 
          )
ssTest2 = @instantiateModel(SSTest2, logCode=false)  
simulate!(ssTest2, stopTime=1.0, log=true, logStates=true)
printResultInfo(ssTest2)

plot(ssTest2, [("submodel.ss.x", "submodel.x1", "submodel.x2", "submodel.x3", "submodel.x4", "ss.x" ),
               ("submodel.ss.u", "ss.u", "y1", "y2"),
               ("submodel.ss.der(x)", "submodel.der(x1)", "submodel.der(x2)", "submodel.der(x3)", "submodel.der(x4)", "ss.der(x)")], figure=3)

simulate!(ssTest2, stopTime=1.0, log=false, logStates=true,
          merge = Map(submodel = Map(x4 = [0.85]), 
                      ss = LinearStateSpace(A=[-1/T   0.0   0.0;
                                               0.0   -1/T   0.0;
                                               0.0    0.0  -1/T],
                                               B=[1.0/T;
                                                  1.0/T;
                                                  1.0/T;;],
                                               C=[0.5 0.5 0.5;],
                                               x_init=[0.35,0.45,0.55]))
         ) 
printResultInfo(ssTest2)

println("\n... Check functions for parameters and signals")
showParameters(         ssTest2)
showEvaluatedParameters(ssTest2)
@test hasParameter(         ssTest2, "submodel.ss.A")
@test getParameter(         ssTest2, "ss.x_init") == [0.35,0.45,0.55]
@test getEvaluatedParameter(ssTest2, "ss.x_init") == [0.35,0.45,0.55]
@test getLastValue(         ssTest2, "ss.x_init") == [0.35,0.45,0.55]

plot(ssTest2, [("submodel.ss.x", "submodel.x1", "submodel.x2", "submodel.x3", "submodel.x4", "ss.x" ),
               ("submodel.ss.u", "ss.u", "y1", "y2"),
               ("submodel.ss.der(x)", "submodel.der(x1)", "submodel.der(x2)", "submodel.der(x3)", "submodel.der(x4)", "ss.der(x)")], maxLegend=12, figure=4)

end